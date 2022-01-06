use std::f64;
use indicatif::ProgressBar;

use mymatrix;
use needleman::Alignment;
use needleman::Direction;
use needleman::Direction::{Diag, Done, Left, Up};
use needleman::Scores;

/// Aligns two sequences using the Needleman Wunsch global alignment with simple gap scoring
pub fn smith_waterman_no_diag(seq1: &Vec<char>, seq2: &Vec<char>, scores: &Scores, min_diag_distance: i32) -> Alignment {
    let seq1_limit = seq1.len() + 1;
    let seq2_limit = seq2.len() + 1;

    let mut mtx = mymatrix::MyMatrix::new(seq1_limit, seq2_limit, 0.0);
    println!("Created an alignment matrix of size [{},{}]", mtx.rows(), mtx.cols());
    let mut trc: mymatrix::MyMatrix<Direction> = mymatrix::MyMatrix::new(seq1_limit, seq2_limit, Done);
    smith_waterman_no_diag_borrow(seq1, seq2, &mut mtx, &mut trc, scores, min_diag_distance)
}

pub fn smith_waterman_no_diag_borrow(seq1: &Vec<char>,
                                     seq2: &Vec<char>,
                                     mtx: &mut mymatrix::MyMatrix<f64>,
                                     trc: &mut mymatrix::MyMatrix<Direction>,
                                     scores: &Scores,
                                     min_diag_distance: i32) -> Alignment {
    let seq1_limit = seq1.len() + 1;
    let seq2_limit = seq2.len() + 1;

    // first square
    mtx.set(0, 0, 0.0);

    // initialize the top row and first column
    for n in 1..seq1_limit {
        mtx.set(n, 0, 0.0);
        trc.set(n, 0, Up);
    }
    for n in 1..seq2_limit {
        mtx.set(0, n, 0.0);
        trc.set(0, n, Left);
    }

    let mut top_score = 0.0;
    let mut topx = 0;
    let mut topy = 0;

    println!("Aligning (status by rows)...");
    let bar = indicatif::ProgressBar::new(seq1_limit as u64);

    // fill in the matrix
    for ix in 1..seq1_limit {
        for iy in 1..seq2_limit {
            let score = Scores::scoring_function(seq1[ix - 1], seq2[iy - 1], scores);
            //if ((ix == iy || ((iy % seq1.len()) == ix) || ((ix % seq2.len()) == iy)) && ix % 100 == 0) {
            //    println!("NOGO is from {},{}", ix, iy);
            //}
            let up_t = (mtx.get(ix - 1, iy) + scores.gap_ext, Up);
            let left_t = (mtx.get(ix, iy - 1) + scores.gap_ext, Left);
            let diag_t = (mtx.get(ix - 1, iy - 1) + score, Diag);

            let mut max = max2(max2(max2(up_t, left_t), diag_t), (0.0, Diag));
            if max.0 > top_score {
                top_score = max.0;
                topx = ix;
                topy = iy;
            }
            if (ix as i32 - iy as i32 ).abs() < min_diag_distance || ((iy as i32 % seq1.len()as i32 ) - ix as i32).abs() < min_diag_distance || ((ix as i32 % seq2.len() as i32) - iy as i32).abs() < min_diag_distance {
                max = (0.0, max.1);
            }

            mtx.set(ix, iy, max.0);
            trc.set(ix, iy, max.1);
        }
        bar.inc(1);
    }
    bar.finish();
    //trc.print_matrix(8);
    //println!("max isnow from {},{}", topx, topy);
    traceback(seq1, seq2, trc, mtx, top_score, topx, topy)
}

#[inline]
fn max2(x: (f64, Direction), y: (f64, Direction)) -> (f64, Direction) {
    if x.0 > y.0 { x } else { y }
}

// a convenience struct for matching movement tuples
#[allow(dead_code)]
struct TruthSet {
    same_row: bool,
    same_col: bool,
    on_diag: bool,
}

/// traceback a matrix into an alignment struct
pub fn traceback(seq1: &Vec<char>,
                 seq2: &Vec<char>,
                 trc: &mymatrix::MyMatrix<Direction>,
                 mtx: &mymatrix::MyMatrix<f64>,
                 top_score: f64,
                 topx: usize,
                 topy: usize) -> Alignment {
    assert_eq!(seq1.len(), trc.rows() - 1, "The matrix doesn't have the right number of rows; rows: {}, expected: {}", seq1.len(), trc.rows() - 1);
    assert_eq!(seq2.len(), trc.cols() - 1, "The matrix doesn't have the right number of columns; columns: {}, expected: {}", seq2.len(), trc.cols() - 1);

    let mut alignment1: Vec<char> = Vec::new();
    let mut alignment2: Vec<char> = Vec::new();

    let mut row_index = topx;
    let mut column_index = topy;

    //println!("Tops: {},{} score {}", topx, topy, top_score);

    let gap = '-';
    //trc.print_matrix(4);

    let mut current_pointer = trc.get(row_index, column_index);
    let mut current_score = mtx.get(row_index, column_index);
    loop {
        match (current_pointer, current_score) {
            (_, x) if x <= 0.0 => {
                //println!("DONE 0: {},{} score {}",row_index,column_index,current_score);
                break;
            }
            (Up, _) => {
                alignment1.push(seq1[row_index - 1]);
                alignment2.push(gap);
                //println!("UP: Moving from {},{} to {},{} score {}", row_index, column_index, row_index - 1, column_index, current_score);
                row_index -= 1;
            }
            (Left, _) => {
                //println!("LEFT: Moving from {},{} to {},{} score {}", row_index, column_index, row_index, column_index - 1, current_score);

                alignment1.push(gap);
                alignment2.push(seq2[column_index - 1]);
                column_index -= 1;
            }
            (Diag, _) => {
                alignment1.push(seq1[row_index - 1]);
                alignment2.push(seq2[column_index - 1]);
                //println!("DIAG: Moving from {},{} to {},{} score {}", row_index, column_index, row_index - 1, column_index - 1, current_score);
                row_index -= 1;
                column_index -= 1;
            }
            (Done, _) => {
                //println!("DONE: {},{} score {}",row_index,column_index,current_score);
                break;
            }
        }
        //println!("looooop");
        current_pointer = trc.get(row_index, column_index);
        current_score = mtx.get(row_index, column_index);
    }

    alignment1.reverse();
    alignment2.reverse();
    println!("Alignment lengths of {} and {}", alignment1.len(), alignment2.len());
    return Alignment {
        seq_one: seq1.to_vec(),
        seq_two: seq2.to_vec(),
        score: top_score,
        start_x: row_index,
        start_y: column_index,
        end_x: topx,
        end_y: topy,
        seq_one_aligned: alignment1,
        seq_two_aligned: alignment2,
    };
}

#[cfg(test)]
mod tests {
    // Note this useful idiom: importing names from outer (for mod tests) scope.
    use super::*;

    #[test]
    fn test_basic_alignment() {
        let scores = Scores::default_scores();
        let alignment = smith_waterman_no_diag(&vec!['A', 'A', 'A', 'A', 'A', 'A'], &vec!['A', 'A', 'A', 'A', 'A', 'A'], &scores, 1);

        let str1: String = alignment.seq_one.into_iter().collect();
        let str2: String = alignment.seq_two.into_iter().collect();
        let str1align: String = alignment.seq_one_aligned.into_iter().collect();
        let str2align: String = alignment.seq_two_aligned.into_iter().collect();
        println!("{},{}", str1align, str2align);

        assert_eq!(alignment.score, 5.0 * scores.match_score);
        assert_eq!(str1, "AAAAAA");
        assert_eq!(str2, "AAAAAA");
        assert_eq!(str1align, "AAAAA");
        assert_eq!(str2align, "AAAAA");
    }

    #[test]
    fn test_basic_alignment_forced() {
        let scores = Scores::default_scores();
        let alignment = smith_waterman_no_diag(&vec!['A', 'C', 'G', 'T', 'A', 'C'], &vec!['A', 'C', 'G', 'T', 'A', 'C'], &scores, 1);

        let str1: String = alignment.seq_one.into_iter().collect();
        let str2: String = alignment.seq_two.into_iter().collect();
        let str1align: String = alignment.seq_one_aligned.into_iter().collect();
        let str2align: String = alignment.seq_two_aligned.into_iter().collect();
        println!("{},{}", str1align, str2align);

        assert_eq!(alignment.score, 2.0 * scores.match_score);
        assert_eq!(str1, "ACGTAC");
        assert_eq!(str2, "ACGTAC");
        assert_eq!(str1align, "AC");
        assert_eq!(str2align, "AC");
    }
}
