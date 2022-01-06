use gcd::Gcd;

/// Transforms a matrix into row echelon form by Gaussian elimination and returns its rank.
pub fn row_reduce(mat: &mut Vec<Vec<i32>>) -> usize {
    let row_n = mat.len();
    let col_n = mat[0].len();

    let mut row_i = 0;
    let mut col_i = 0;

    while row_i < row_n && col_i < col_n {
        // Locate a pivot with the least absolute value.
        let pivot = mat
            .iter()
            .map(|row| row[col_i])
            .enumerate()
            .skip(row_i)
            .filter(|(_, n)| *n != 0)
            .min_by_key(|(_, n)| n.unsigned_abs());

        let (prev_pivot_i, pivot) = match pivot {
            Some(p) => p,
            // All zero.
            None => {
                col_i += 1;
                continue;
            }
        };

        mat.swap(row_i, prev_pivot_i);
        let (pivot_row, rest) = mat[row_i..].split_first_mut().unwrap();

        for cur_row in rest {
            let cur = cur_row[col_i];
            if cur == 0 {
                continue;
            }

            let (cur_m, pivot_m) = lcm_div(cur, pivot);

            row_mul_sub(cur_row, cur_m, pivot_row, pivot_m);
        }

        row_i += 1;
        col_i += 1;
    }
    row_i
}

/// Returns the positive GCD of two signed integers.
fn gcd(u: i32, v: i32) -> i32 {
    u.unsigned_abs().gcd(v.unsigned_abs()) as i32
}

/// Returns the multipliers that, when multiplied by the
/// corresponding integer, give the LCM of two integers.
fn lcm_div(a: i32, b: i32) -> (i32, i32) {
    let gcd = gcd(a, b);
    (b / gcd, a / gcd)
}

fn row_mul_sub(dest: &mut [i32], dest_m: i32, src: &[i32], src_m: i32) {
    dest.iter_mut().zip(src.iter()).for_each(|(dest, src)| {
        *dest = *dest * dest_m - *src * src_m;
    })
}

/// Solves the null space of a matrix in row echelon form and returns one possible basis.
pub fn solve(mat: Vec<Vec<i32>>, rank: usize) -> Vec<Vec<i32>> {
    let col_n = mat[0].len();
    match col_n - rank {
        0 => vec![],
        1 => vec![solve_unique(mat, rank)],
        _ => solve_multiple(mat, rank),
    }
}

fn solve_unique(mat: Vec<Vec<i32>>, rank: usize) -> Vec<i32> {
    let mut sol = vec![0; rank + 1];
    sol[rank] = 1;

    for i in (0..rank).rev() {
        let sum: i32 = mat[i][i + 1..]
            .iter()
            .zip(sol[i + 1..].iter())
            .map(|(&a, &b)| a * b)
            .sum();
        let pivot = mat[i][i];

        let (mut sum_m, mut pivot_m) = lcm_div(sum, pivot);

        if sum_m < 0 {
            sum_m = -sum_m;
        } else {
            pivot_m = -pivot_m;
        }

        sol[i + 1..].iter_mut().for_each(|n| *n *= sum_m);
        sol[i] = pivot_m;
    }

    normalize_sol(&mut sol);
    sol
}

/// Normalizes a solution, so that the number of positive coefficients
/// are greater than that of negative ones or, if the numbers are equal,
/// the first nonzero coefficient is positive.
fn normalize_sol(sol: &mut [i32]) {
    let mut neg: i32 = sol.iter().map(|n| n.signum()).sum();
    if neg == 0 {
        neg = sol.iter().copied().find(|n| *n != 0).unwrap();
    }
    if neg < 0 {
        sol.iter_mut().for_each(|n| *n = -*n);
    }
}

fn solve_multiple(mut mat: Vec<Vec<i32>>, rank: usize) -> Vec<Vec<i32>> {
    let col_n = mat[0].len();
    let nullity = col_n - rank;

    row_augment_identity(&mut mat, rank);
    col_reduce_upper(&mut mat, rank);

    // Transpose and move the basis to the top, without extra allocation.
    for i in 0..col_n.min(nullity) {
        // Move the outermost row and column at a time, so that
        // no vector is overwritten before they are moved.
        for j in i..col_n {
            mat[i][j] = mat[rank + j][rank + i];
        }
        for j in i + 1..nullity {
            mat[j][i] = mat[rank + i][rank + j];
        }
    }

    mat.truncate(nullity);
    for sol in &mut mat {
        normalize_sol(sol);
        simplify_sol(sol);
    }
    mat
}

/// Divides the GCD from a solution.
fn simplify_sol(sol: &mut [i32]) {
    let gcd = sol.iter().copied().reduce(gcd).unwrap();
    if gcd != 1 {
        sol.iter_mut().for_each(|n| *n /= gcd);
    }
}

fn row_augment_identity(mat: &mut Vec<Vec<i32>>, rank: usize) {
    let row_n = mat.len();
    let col_n = mat[0].len();

    mat.extend((row_n - rank..col_n).map(|_| vec![0; col_n]));

    mat[rank..]
        .iter_mut()
        .take(col_n)
        .enumerate()
        .for_each(|(i, row)| row[i] = 1);
}

fn col_reduce_upper(mat: &mut Vec<Vec<i32>>, upper_rank: usize) {
    let col_n = mat[0].len();

    for i in 0..upper_rank {
        let (prev_pivot_i, pivot) = mat[i]
            .iter()
            .copied()
            .enumerate()
            .skip(i)
            .find(|(_, n)| *n != 0)
            .unwrap(); // There must be a pivot.

        col_swap(mat, i, prev_pivot_i);
        let pivot_i = i;

        for cur_i in pivot_i + 1..col_n {
            let cur = mat[pivot_i][cur_i];
            if cur == 0 {
                continue;
            }

            let (cur_m, pivot_m) = lcm_div(cur, pivot);

            col_mul_sub(mat, cur_i, cur_m, pivot_i, pivot_m);
        }
    }
}

fn col_swap(mat: &mut Vec<Vec<i32>>, a: usize, b: usize) {
    if a != b {
        for row in mat {
            row.swap(a, b);
        }
    }
}

fn col_mul_sub(mat: &mut Vec<Vec<i32>>, dest_i: usize, dest_m: i32, src_i: usize, src_m: i32) {
    for row in mat {
        row[dest_i] = row[dest_i] * dest_m - row[src_i] * src_m;
    }
}
