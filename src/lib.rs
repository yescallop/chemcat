pub fn row_reduce(mat: &mut Vec<Vec<i32>>) {
    let w = mat[0].len();
    let h = mat.len();

    for i in 0..w.min(h) {
        // Locate a pivot with the least absolute value.
        let mut min = (0, u32::MAX);
        for j in i..h {
            let cur = mat[j][i];
            if cur == 0 {
                continue;
            }
            let cur_abs = cur.abs() as u32;
            if cur_abs < min.1 {
                min = (j, cur_abs);
            }
        }

        if min.1 == u32::MAX {
            break;
        }

        mat.swap(i, min.0);
        let pivot = mat[i][i];
        let (pivot_row, rest) = mat[i..].split_first_mut().unwrap();

        for cur_row in rest {
            let cur = cur_row[i];
            if cur == 0 {
                continue;
            }

            let gcd = gcd(cur, pivot);
            let cur_m = pivot / gcd;
            let pivot_m = cur / gcd;

            mul_sub(cur_row, cur_m, pivot_row, pivot_m);
        }
    }
}

fn mul_sub(dest: &mut [i32], dest_m: i32, src: &[i32], src_m: i32) {
    dest.iter_mut().zip(src.iter()).for_each(|(dest, src)| {
        *dest = *dest * dest_m - *src * src_m;
    })
}

fn gcd(mut a: i32, mut b: i32) -> i32 {
    while b != 0 {
        let t = b;
        b = a % b;
        a = t;
    }
    a
}

pub fn solve(mut mat: Vec<Vec<i32>>) -> Solution {
    let zero_rows = mat
        .iter()
        .rev()
        .take_while(|row| row.iter().all(|cur| *cur == 0))
        .count();
    let nonzero_rows = mat.len() - zero_rows;
    let unknowns = mat[0].len();
    match unknowns - nonzero_rows {
        0 => Solution::None,
        1 => {
            mat.truncate(nonzero_rows);
            Solution::Unique(solve_unique(mat))
        }
        n => Solution::Infinite(n as u32),
    }
}

fn solve_unique(mat: Vec<Vec<i32>>) -> Vec<u32> {
    let h = mat.len();
    let mut sol = vec![0; h + 1];
    sol[h] = 1;

    for i in (0..h).rev() {
        let sum: i32 = mat[i][i + 1..]
            .iter()
            .zip(sol[i + 1..].iter())
            .map(|(&a, &b)| a * b as i32)
            .sum();
        let pivot = mat[i][i];

        let gcd = gcd(sum, pivot);
        let sum_m = (pivot / gcd).abs() as u32;
        let pivot_m = (sum / gcd).abs() as u32;

        sol[i + 1..].iter_mut().for_each(|n| *n *= sum_m);
        sol[i] = pivot_m;
    }

    sol
}

pub enum Solution {
    None,
    Unique(Vec<u32>),
    Infinite(u32),
}
