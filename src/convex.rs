use mymatrix;
use needleman::Scores;
use needleman::Alignment;


/// Aligns two sequences using the Needleman Wunsch global alignment with simple gap scoring
pub fn convex(seq1: &Vec<char>, seq2: &Vec<char>, scores: &Scores) -> Alignment {
    let seq1_limit = seq1.len() + 1;
    let seq2_limit = seq2.len() + 1;

    let mut mtx = mymatrix::MyMatrix::new(seq1_limit, seq2_limit, 0.0);
    ////println!("made matrix of {} {}", mtx.rows(), mtx.cols());
    let mut trc: mymatrix::MyMatrix<i32> = mymatrix::MyMatrix::new(seq1_limit, seq2_limit, 0);
    convex_alignment(seq1, seq2, &mut mtx, &mut trc, scores)
}

pub fn convex_alignment(seq1: &Vec<char>,
                        seq2: &Vec<char>,
                        mtx: &mut mymatrix::MyMatrix<f64>,
                        trc: &mut mymatrix::MyMatrix<i32>,
                        scores: &Scores) -> Alignment {

    let seq1_limit = seq1.len() + 1;
    let seq2_limit = seq2.len() + 1;
    assert!(mtx.rows() == seq1_limit);
    assert!(mtx.cols() == seq2_limit);
    assert!(trc.rows() == seq1_limit);
    assert!(trc.cols() == seq2_limit);

    // TRC is stored as offsets
    // Negative traceback values mean to the left; positive traceback means up,
    // and zero means diagonal (one base)

    // first square
    mtx.set(0, 0, 0.0);

    let scoring_function = |i: usize| -> f64 { -10.0 - (1.0 * ((i + 1) as f64).log2()) };

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
            let up = mymatrix::maximize_over_column(&mtx, iy, ix, &scoring_function);
            let left = mymatrix::maximize_over_row(&mtx, ix, iy, &scoring_function);
            let diag = mtx.get(ix - 1, iy - 1) + score;

            if up.1 > left.1 {
                if diag < up.1 {
                    mtx.set(ix, iy, up.1);
                    trc.set(ix, iy, (ix - up.0) as i32);
                } else {
                    mtx.set(ix, iy, diag);
                    trc.set(ix, iy, 0);
                }
            } else {
                if diag < left.1 {
                    mtx.set(ix, iy, left.1);
                    trc.set(ix, iy, -1 * (iy as i32 - left.0 as i32));
                } else {
                    mtx.set(ix, iy, diag);
                    trc.set(ix, iy, 0);
                }
            };
        }
    }

    let no_cost = |_i: usize| -> f64 {0 as f64};
    //println!("row_index={},rows={} -- {},{}",mtx.rows(),mtx.cols(),seq1_limit,seq2_limit);
    let start_row = mymatrix::maximize_over_column(&mtx, seq2_limit - 1, seq1_limit - 1, &no_cost);
    let start_column = mymatrix::maximize_over_row(&mtx, seq1_limit - 1, seq2_limit - 1, &no_cost);
    let lower_right = mtx.get(seq1_limit - 1, seq2_limit - 1);

    if lower_right > start_column.1 {
        if lower_right > start_row.1 {
            traceback(seq1, seq2, seq1_limit - 1, seq2_limit - 1, trc, mtx.get(seq1_limit - 1, seq2_limit - 1))
        } else {
            traceback(seq1, seq2, start_row.0, seq2_limit - 1, trc, mtx.get(seq1_limit - 1, seq2_limit - 1))
        }
    } else {
        if start_column.1 > start_row.1 {
            traceback(seq1, seq2, seq1_limit - 1, start_column.0, trc, mtx.get(seq1_limit - 1, seq2_limit - 1))
        } else {
            traceback(seq1, seq2, start_row.0, seq2_limit - 1, trc, mtx.get(seq1_limit - 1, seq2_limit - 1))
        }
    }
}

/// traceback a matrix into an alignment struct
pub fn traceback(seq1: &Vec<char>, seq2: &Vec<char>, start_row: usize, start_column: usize, trc: &mymatrix::MyMatrix<i32>, score: f64) -> Alignment {
    assert_eq!(seq1.len(), trc.rows() - 1, "The matrix doesn't have the right number of rows; rows: {}, expected: {}", seq1.len(), trc.rows() - 1);
    assert_eq!(seq2.len(), trc.cols() - 1, "The matrix doesn't have the right number of columns; columns: {}, expected: {}", seq2.len(), trc.cols() - 1);

    let mut alignment1 = Vec::new();
    let mut alignment2 = Vec::new();

    let mut row_index = start_row as u32;
    let mut column_index = start_column as u32;

    // if we're off the final score
    if row_index + 1 < trc.rows() as u32 {
        alignment2.append(&mut gap_of_length((trc.rows() - 1) - row_index as usize));
        let alignment1_reversed = &mut seq1[row_index as usize ..((trc.rows() - 1) as usize)].to_vec();
        alignment1_reversed.reverse();
        alignment1.append(alignment1_reversed);
    }
    if column_index + 1 < trc.cols() as u32 {
        let alignment2_reversed = &mut seq2[column_index as usize ..((trc.cols() - 1) as usize)].to_vec();
        alignment2_reversed.reverse();
        alignment2.append(alignment2_reversed);
        alignment1.append(&mut gap_of_length((trc.cols() - 1) - column_index as usize));
    }

    let mut current_pointer = trc.get(row_index as usize, column_index as usize);

    loop {
        match current_pointer {
            _ if row_index == 0 && column_index == 0 => {
                break;
            }
            0 => {
                assert!(row_index > 0);
                assert!(column_index > 0);
                alignment1.push(seq1[(row_index - 1) as usize ]);
                alignment2.push(seq2[(column_index -1) as usize ]);
                row_index -= 1;
                column_index -= 1;
            }
            _x if _x < 0 => {
                let offset = (-1 * _x) as usize;
                let move_to_column = column_index as usize - offset;

                assert!(((-1 * _x) as u32) <= column_index);
                alignment1.append(&mut gap_of_length(offset));
                let alignment2_reversed = &mut seq2[move_to_column..(column_index as usize)].to_vec();
                alignment2_reversed.reverse();
                alignment2.append(alignment2_reversed);
                column_index = move_to_column as u32;
            }
            _x if _x > 0 => {
                // _x is the offset, determine the row we're moving to
                let move_to_row = row_index - (_x as u32);
                assert!((_x as u32) <= row_index);

                let alignment1_reversed = &mut seq1[(move_to_row as usize)..(row_index as usize)].to_vec();
                alignment1_reversed.reverse();
                alignment1.append(alignment1_reversed);
                alignment2.append(&mut gap_of_length(_x as usize));
                row_index = move_to_row as u32;
            }
            _ => unreachable!()
        }
        current_pointer = trc.get(row_index as usize, column_index as usize);
    }
    alignment1.reverse();
    alignment2.reverse();

    return Alignment {
        seq_one: seq1.to_vec(),
        seq_two: seq2.to_vec(),
        score: score,
        seq_one_aligned: alignment1,
        seq_two_aligned: alignment2,
    };
}

fn gap_of_length(x: usize) -> Vec<char> {
    (0..x).map(|_| '-').collect::<Vec<char>>()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_basic_alignment() {
        let scores = Scores::default_scores();
        let mut mtx = mymatrix::MyMatrix::new(4, 4, 0.0);
        let mut trc: mymatrix::MyMatrix<i32> = mymatrix::MyMatrix::new(4, 4, 0);

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


    #[test]
    fn test_bigger_gap_alignment() {
        let scores = Scores::default_scores();
        let seq1 = "AAATTTTTTTTTTTTTTTTTTTTTTTTAAA";
        let seq2 = "AAAAAA";
        let mut mtx = mymatrix::MyMatrix::new(seq1.len() + 1, seq2.len() + 1, 0.0);
        let mut trc: mymatrix::MyMatrix<i32> = mymatrix::MyMatrix::new(seq1.len() + 1, seq2.len() + 1, 0);

        let alignment = convex_alignment(&seq1.chars().collect::<Vec<char>>(),
                                         &seq2.chars().collect::<Vec<char>>(),
                                         &mut mtx,
                                         &mut trc,
                                         &scores);

        let str1: String = alignment.seq_one.into_iter().collect();
        let str2: String = alignment.seq_two.into_iter().collect();
        let str1align: String = alignment.seq_one_aligned.into_iter().collect();
        let str2align: String = alignment.seq_two_aligned.into_iter().collect();

        assert_eq!(str1align, "AAATTTTTTTTTTTTTTTTTTTTTTTTAAA");
        assert_eq!(str2align, "AAA------------------------AAA");
    }

    #[test]
    fn test_bigger_gap_alignment_reversed() {
        let scores = Scores::default_scores();
        let seq2 = "AAATTTTTTTTTTTTTTTTTTTTTTTTAAA";
        let seq1 = "AAAAAA";
        let mut mtx = mymatrix::MyMatrix::new(seq1.len() + 1, seq2.len() + 1, 0.0);
        let mut trc: mymatrix::MyMatrix<i32> = mymatrix::MyMatrix::new(seq1.len() + 1, seq2.len() + 1, 0);

        let alignment = convex_alignment(&seq1.chars().collect::<Vec<char>>(),
                                         &seq2.chars().collect::<Vec<char>>(),
                                         &mut mtx,
                                         &mut trc,
                                         &scores);

        let str1: String = alignment.seq_one.into_iter().collect();
        let str2: String = alignment.seq_two.into_iter().collect();
        let str1align: String = alignment.seq_one_aligned.into_iter().collect();
        let str2align: String = alignment.seq_two_aligned.into_iter().collect();

        assert_eq!(str2align, "AAATTTTTTTTTTTTTTTTTTTTTTTTAAA");
        assert_eq!(str1align, "AAA------------------------AAA");
    }

    #[test]
    fn test_end_with_dels_seq1() {
        let scores = Scores::default_scores();
        let seq2 = "AATAAAGGTGGG";
        let seq1 = "AATAAA";
        let mut mtx = mymatrix::MyMatrix::new(seq1.len() + 1, seq2.len() + 1, 0.0);
        let mut trc: mymatrix::MyMatrix<i32> = mymatrix::MyMatrix::new(seq1.len() + 1, seq2.len() + 1, 0);

        let alignment = convex_alignment(&seq1.chars().collect::<Vec<char>>(),
                                         &seq2.chars().collect::<Vec<char>>(),
                                         &mut mtx,
                                         &mut trc,
                                         &scores);

        let str1: String = alignment.seq_one.into_iter().collect();
        let str2: String = alignment.seq_two.into_iter().collect();
        let str1align: String = alignment.seq_one_aligned.into_iter().collect();
        let str2align: String = alignment.seq_two_aligned.into_iter().collect();
        //println!("1 = {}",str1align);
        //println!("2 = {}",str2align);
        assert_eq!(str2align, "AATAAAGGTGGG");
        assert_eq!(str1align, "AATAAA------");
    }

    #[test]
    fn test_end_with_dels_seq2() {
        let scores = Scores::default_scores();
        let seq2 = "AAAAAAGGGGGGGGGGGGGGGGGGGGGGTTT";
        let seq1 = "GGGTTT";
        let mut mtx = mymatrix::MyMatrix::new(seq1.len() + 1, seq2.len() + 1, 0.0);
        let mut trc: mymatrix::MyMatrix<i32> = mymatrix::MyMatrix::new(seq1.len() + 1, seq2.len() + 1, 0);

        let alignment = convex_alignment(&seq1.chars().collect::<Vec<char>>(),
                                         &seq2.chars().collect::<Vec<char>>(),
                                         &mut mtx,
                                         &mut trc,
                                         &scores);

        let str1: String = alignment.seq_one.into_iter().collect();
        let str2: String = alignment.seq_two.into_iter().collect();
        let str1align: String = alignment.seq_one_aligned.into_iter().collect();
        let str2align: String = alignment.seq_two_aligned.into_iter().collect();

        assert_eq!(str2align, "AAAAAAGGGGGGGGGGGGGGGGGGGGGGTTT");
        assert_eq!(str1align, "-------------------------GGGTTT");
    }
}