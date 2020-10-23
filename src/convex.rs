use mymatrix;
use needleman::Scores;
use needleman::Alignment;


pub fn convex_alignment(seq1: &Vec<char>,
                        seq2: &Vec<char>,
                        mtx: &mut mymatrix::MyMatrix<f64>,
                        trc: &mut mymatrix::MyMatrix<i32>,
                        scores: &Scores) -> Alignment {
    let seq1_limit = seq1.len() + 1;
    let seq2_limit = seq2.len() + 1;

    // Negative traceback values mean to the left; positive traceback means up,
    // and zero means diagonal (one base)

    // first square
    mtx.set(0, 0, 0.0);

    let scoring_function = |i: usize| -> f64 { -10.0 + (-1.0 * ((i + 1) as f64).log2()) };

    // initialize the first column and top row
    for n in 1..seq1_limit {
        mtx.set(n, 0, scoring_function(n));
        trc.set(n, 0, n as i32);
    }
    for n in 1..seq2_limit {
        mtx.set(0, n, scoring_function(n));
        trc.set(0, n, -1 * (n as i32));
    }

    // fill in the matrix
    for ix in 1..seq1_limit {
        for iy in 1..seq2_limit {
            let score = Scores::scoring_function(seq1[ix - 1], seq2[iy - 1], scores);
            let up = mymatrix::maximize_over_column(&mtx, ix, iy, &scoring_function);
            let left = mymatrix::maximize_over_row(&mtx, ix, iy, &scoring_function);
            let diag = mtx.get(ix - 1, iy - 1) + score;

            if up.1 > left.1 {
                if diag < up.1 {
                    mtx.set(ix, iy, up.1);
                    trc.set(ix, iy, up.0 as i32);
                } else {
                    mtx.set(ix, iy, diag);
                    trc.set(ix, iy, 0);
                }
            } else {
                if diag < left.1 {
                    mtx.set(ix, iy, left.1);
                    trc.set(ix, iy, -1 * (left.0 as i32));
                } else {
                    mtx.set(ix, iy, diag);
                    trc.set(ix, iy, 0);
                }
            };
        }
    }

    traceback(seq1, seq2, trc, mtx.get(seq1_limit - 1, seq2_limit - 1))
}

/// traceback a matrix into an alignment struct
pub fn traceback(seq1: &Vec<char>, seq2: &Vec<char>, trc: &mymatrix::MyMatrix<i32>, score: f64) -> Alignment {
    assert_eq!(seq1.len(), trc.rows() - 1, "The matrix doesn't have the right number of rows; rows: {}, expected: {}", seq1.len(), trc.rows() - 1);
    assert_eq!(seq2.len(), trc.cols() - 1, "The matrix doesn't have the right number of columns; columns: {}, expected: {}", seq2.len(), trc.cols() - 1);

    let mut alignment1 = Vec::new();
    let mut alignment2 = Vec::new();

    let mut row_index = seq1.len() as i32;
    let mut column_index = seq2.len() as i32;

    let gap = '-';
    //trc.print_matrix();

    let mut current_pointer = trc.get(row_index as usize, column_index as usize);

    loop {
        // neg = left, positive = up
        match current_pointer {
            0 => {
                alignment1.push(seq1[row_index as usize - 1]);
                alignment2.push(seq2[column_index as usize - 1]);
                // println!("DIAG: Moving from {},{} to {},{}", row_index, column_index, row_index - 1, column_index - 1);
                row_index -= 1;
                column_index -= 1;
            }
            _x if _x < 0 => {
                alignment1.append(&mut gap_of_length((row_index - (-1 * _x)) as usize));
                alignment2.append(&mut seq2[(_x as usize)..((column_index - 1) as usize)].to_vec());
                // println!("LEFT: Moving from {},{} to {},{}", row_index, column_index, row_index, column_index - 1);
                row_index += _x;
            }
            _x if _x > 0 => {
                alignment1.append(&mut seq1[(_x as usize)..(row_index as usize - 1)].to_vec());
                alignment2.append(&mut gap_of_length((column_index - _x) as usize));
                // println!("LEFT: Moving from {},{} to {},{}", row_index, column_index, row_index - 1, column_index);
                column_index -= _x;
            }
            _ => unreachable!()
        }
        current_pointer = trc.get(row_index as usize, column_index as usize);
    }

    alignment1.reverse();
    alignment2.reverse();
    // println!("{},{}",alignment1.len(),alignment2.len());
    return Alignment {
        seq_one: seq1.to_vec(),
        seq_two: seq2.to_vec(),
        score: score,
        seq_one_aligned: alignment1,
        seq_two_aligned: alignment2,
    };
}

fn gap_of_length(x: usize) -> Vec<char> {
    (0..10).map(|_| '-').collect::<Vec<char>>()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_basic_alignment() {
        let scores = Scores::default_scores();
        let mut mtx = mymatrix::MyMatrix::new(3, 3, 0.0);
        let mut trc: mymatrix::MyMatrix<i32> = mymatrix::MyMatrix::new(3, 3, 0);

        let alignment = convex_alignment(&vec!['A', 'A', 'A'],
                                         &vec!['A', 'A', 'A'],
                                         &mut mtx,
                                         &mut trc,
                                         &scores);

        let str1: String = alignment.seq_one.into_iter().collect();
        let str2: String = alignment.seq_two.into_iter().collect();
        let str1align: String = alignment.seq_one_aligned.into_iter().collect();
        let str2align: String = alignment.seq_two_aligned.into_iter().collect();

        assert_eq!(str1, "AAA");
        assert_eq!(str2, "AAA");
        assert_eq!(str1align, "AAA");
        assert_eq!(str2align, "AAA");
        assert_eq!(alignment.score, 3.0 * scores.match_score);
    }
}