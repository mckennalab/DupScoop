use mymatrix;
use needleman::Direction::{Diag, Done, Left, Up};
use std::f64;
use std::fmt;
use needleman::{Direction, Scores, Alignment};

/// Aligns two sequences using the Needleman Wunsch global alignment with simple gap scoring
pub fn affine_align(seq1: &Vec<char>, seq2: &Vec<char>, scores: &Scores) -> Alignment {
    let seq1_limit = seq1.len() + 1;
    let seq2_limit = seq2.len() + 1;

    let mut match_matrix = mymatrix::MyMatrix::new(seq1_limit, seq2_limit, 0.0);
    let mut ins_matrix = mymatrix::MyMatrix::new(seq1_limit, seq2_limit, 0.0);
    let mut del_matrix = mymatrix::MyMatrix::new(seq1_limit, seq2_limit, 0.0);

    println!("made matrix of {} {}", match_matrix.rows(), match_matrix.cols());
    let mut match_trc: mymatrix::MyMatrix<Direction> = mymatrix::MyMatrix::new(seq1_limit, seq2_limit, Done);
    let mut ins_trc: mymatrix::MyMatrix<Direction> = mymatrix::MyMatrix::new(seq1_limit, seq2_limit, Done);
    let mut del_trc: mymatrix::MyMatrix<Direction> = mymatrix::MyMatrix::new(seq1_limit, seq2_limit, Done);

    affine_borrow(seq1,
                  seq2,
                  &mut match_matrix,
                  &mut ins_matrix,
                  &mut del_matrix,
                  &mut match_trc,
                  &mut ins_trc,
                  &mut del_trc,
                  scores)
}

pub fn affine_borrow(seq1: &Vec<char>,
                     seq2: &Vec<char>,
                     match_matrix: &mut mymatrix::MyMatrix<f64>,
                     ins_matrix: &mut mymatrix::MyMatrix<f64>,
                     del_matrix: &mut mymatrix::MyMatrix<f64>,
                     match_trc: &mut mymatrix::MyMatrix<Direction>,
                     ins_trc: &mut mymatrix::MyMatrix<Direction>,
                     del_trc: &mut mymatrix::MyMatrix<Direction>,
                     scores: &Scores) -> Alignment {
    let seq1_limit = seq1.len() + 1;
    let seq2_limit = seq2.len() + 1;

    // we want a practical MAX, but not at the limit of F64 values
    let practical_max = 10000000000.0;

    // first square
    match_matrix.set(0, 0, practical_max);

    // initialize the top row and first column
    for n in 1..seq1_limit {
        match_matrix.set(n, 0, practical_max);
        ins_matrix.set(n, 0, scores.gap_ext * (n as f64));
        del_matrix.set(n, 0, practical_max);
        trc.set(n, 0, Up);
    }
    for n in 1..seq2_limit {
        match_matrix.set(0, n, practical_max);
        ins_matrix.set(n, 0, practical_max);
        del_matrix.set(n, 0, scores.gap_ext * (n as f64));
        trc.set(0, n, Left);
    }


    // fill in the matrix
    for ix in 1..seq1_limit {
        for iy in 1..seq2_limit {
            let score = Scores::scoring_function(seq1[ix - 1], seq2[iy - 1], scores);

            let match_max = max2(
                    max2((score + ins_matrix.get(ix - 1, iy - 1),Up),
                    (score + del_matrix.get(ix - 1, iy - 1),Left)),
                (score + match_matrix.get(ix - 1, iy - 1),Diag));

            match_matrix.set(ix, iy, match_max.0);
            match_trc.set(ix, iy, match_max.1);

            let ins_max = max2(
                max2((scores.gap_open + scores.gap_ext + match_matrix.get(ix, iy - 1),Diag),
                     (ins_matrix.get(ix, iy - 1),Up)),
                (scores.gap_open + scores.gap_ext + del_matrix.get(ix, iy - 1),Left));

            ins_matrix.set(ix, iy, ins_max.0);
            ins_trc.set(ix, iy, ins_max.1);


            let del_max = max2(
                max2((scores.gap_open + scores.gap_ext + match_matrix.get(ix, iy - 1),Diag),
                     (del_matrix.get(ix, iy - 1),Left)),
                (scores.gap_open + scores.gap_ext + ins_matrix.get(ix, iy - 1),Up));

            del_matrix.set(ix, iy, del_max.0);
            del_trc.set(ix, iy, del_max.1);

            if ix % 1000 == 0 && iy % 1000 == 0 {
                println!("top is now {} from {},{}", match_matrix.0, ix, iy);
            }
        }
    }
    //trc.print_matrix(8);
    traceback(seq1, seq2, trc, mtx.get(seq1_limit - 1, seq2_limit - 1), seq1_limit - 1, seq2_limit - 1)
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
pub fn traceback(seq1: &Vec<char>, seq2: &Vec<char>, trc: &mymatrix::MyMatrix<Direction>, top_score: f64, topx: usize, topy: usize) -> Alignment {
    assert_eq!(seq1.len(), trc.rows() - 1, "The matrix doesn't have the right number of rows; rows: {}, expected: {}", seq1.len(), trc.rows() - 1);
    assert_eq!(seq2.len(), trc.cols() - 1, "The matrix doesn't have the right number of columns; columns: {}, expected: {}", seq2.len(), trc.cols() - 1);

    let mut alignment1 = Vec::new();
    let mut alignment2 = Vec::new();

    let mut row_index = seq1.len(); // topx;
    let mut column_index = seq2.len(); // topx;

    let gap = '-';

    let mut current_pointer = trc.get(row_index, column_index);

    loop {
        match current_pointer {
            Up => {
                alignment1.push(seq1[row_index - 1]);
                alignment2.push(gap);
                // println!("LEFT: Moving from {},{} to {},{}", row_index, column_index, row_index - 1, column_index);
                row_index -= 1;
            }
            Left => {
                alignment1.push(gap);
                alignment2.push(seq2[column_index - 1]);
                // println!("LEFT: Moving from {},{} to {},{}", row_index, column_index, row_index, column_index - 1);
                column_index -= 1;
            }
            Diag => {
                alignment1.push(seq1[row_index - 1]);
                alignment2.push(seq2[column_index - 1]);
                // println!("DIAG: Moving from {},{} to {},{}", row_index, column_index, row_index - 1, column_index - 1);
                row_index -= 1;
                column_index -= 1;
            }
            Done => {
                //println!("{},{}",row_index,column_index);
                break;
            }
        }
        current_pointer = trc.get(row_index, column_index);
    }

    alignment1.reverse();
    alignment2.reverse();
    // println!("{},{}",alignment1.len(),alignment2.len());
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
    fn test_max_2() {
        let up = (30.0, Up);
        let left = (40.0, Left);
        let diag = (50.0, Diag);
        let max = max2(max2(diag, left), up);
        assert_eq!(max.0, 50.0);
        assert_eq!(max.1, Diag);

        let up = (60.0, Up);
        let left = (40.0, Left);
        let diag = (50.0, Diag);
        let max = max2(max2(diag, left), up);
        assert_eq!(max.0, 60.0);
        assert_eq!(max.1, Up);

        let up = (30.0, Up);
        let left = (60.0, Left);
        let diag = (50.0, Diag);
        let max = max2(max2(diag, left), up);
        assert_eq!(max.0, 60.0);
        assert_eq!(max.1, Left);
    }

    #[test]
    fn test_basic_alignment() {
        let scores = Scores::default_scores();
        let alignment = needleman_wunsch(&vec!['A', 'A', 'A'], &vec!['A', 'A', 'A'], &scores);

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
    fn test_basic_alignment_ns() {
        let scores = Scores::default_scores();
        let alignment = needleman_wunsch(&vec!['A', 'N', 'A'], &vec!['A', 'T', 'A'], &scores);

        let str1: String = alignment.seq_one.into_iter().collect();
        let str2: String = alignment.seq_two.into_iter().collect();
        let str1align: String = alignment.seq_one_aligned.into_iter().collect();
        let str2align: String = alignment.seq_two_aligned.into_iter().collect();

        assert_eq!(str1, "ANA");
        assert_eq!(str2, "ATA");
        assert_eq!(str1align, "ANA");
        assert_eq!(str2align, "ATA");
        assert_eq!(alignment.score, 2.0 * scores.match_score);
    }

    #[test]
    fn test_basic_alignment_unequal() {
        let scores = Scores::default_scores();
        let alignment = needleman_wunsch(&vec!['A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A'], &vec!['A', 'A', 'A'], &scores);

        let str1: String = alignment.seq_one.into_iter().collect();
        let str2: String = alignment.seq_two.into_iter().collect();
        let str1align: String = alignment.seq_one_aligned.into_iter().collect();
        let str2align: String = alignment.seq_two_aligned.into_iter().collect();

        assert_eq!(str1, "AAAAAAAAAAAA");
        assert_eq!(str2, "AAA");
        assert_eq!(str1align, "AAAAAAAAAAAA");
        assert_eq!(str2align, "---------AAA");
    }

    #[test]
    fn test_basic_alignment_unequal_second() {
        let scores = Scores::default_scores();
        let alignment = needleman_wunsch(&vec!['A', 'A', 'A'], &vec!['A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A'], &scores);

        let str1: String = alignment.seq_one.into_iter().collect();
        let str2: String = alignment.seq_two.into_iter().collect();
        let str1align: String = alignment.seq_one_aligned.into_iter().collect();
        let str2align: String = alignment.seq_two_aligned.into_iter().collect();

        assert_eq!(str1, "AAA");
        assert_eq!(str2, "AAAAAAAAAAAA");
        assert_eq!(str1align, "---------AAA");
        assert_eq!(str2align, "AAAAAAAAAAAA");
    }

    #[test]
    fn test_mismatch_alignment() {
        let scores = Scores::default_scores();
        let alignment = needleman_wunsch(&vec!['A', 'T', 'A'], &vec!['A', 'A', 'A'], &scores);

        let str1: String = alignment.seq_one.into_iter().collect();
        let str2: String = alignment.seq_two.into_iter().collect();
        let str1align: String = alignment.seq_one_aligned.into_iter().collect();
        let str2align: String = alignment.seq_two_aligned.into_iter().collect();

        assert_eq!(str1, "ATA");
        assert_eq!(str2, "AAA");
        println!("Alignment 1: {} alignment 2: {}", str1align, str2align);
        assert_eq!(str1align, "ATA");
        assert_eq!(str2align, "AAA");
        assert_eq!(alignment.score, 2.0 * scores.match_score + scores.mismatch_score);
    }

    #[test]
    fn test_mismatch_length() {
        let scores = Scores::default_scores();
        let alignment = needleman_wunsch(&vec!['A', 'T', 'T', 'A', 'A'], &vec!['T', 'A', 'A'], &scores);

        let str1: String = alignment.seq_one.into_iter().collect();
        let str2: String = alignment.seq_two.into_iter().collect();
        let str1align: String = alignment.seq_one_aligned.into_iter().collect();
        let str2align: String = alignment.seq_two_aligned.into_iter().collect();

        assert_eq!(str1, "ATTAA");
        assert_eq!(str2, "TAA");
        println!("Alignment 1: {} alignment 2: {}", str1align, str2align);

        assert_eq!(str1align, "ATTAA");
        assert_eq!(str2align, "--TAA");
    }


    #[test]
    fn test_mismatch_length2() {
        let scores = Scores::default_scores();
        let alignment = needleman_wunsch(&vec!['G', 'G', 'G', 'A', 'T', 'T', 'A', 'A'], &vec!['G', 'G', 'T', 'A', 'A'], &scores);

        let str1: String = alignment.seq_one.into_iter().collect();
        let str2: String = alignment.seq_two.into_iter().collect();
        let str1align: String = alignment.seq_one_aligned.into_iter().collect();
        let str2align: String = alignment.seq_two_aligned.into_iter().collect();

        assert_eq!(str1, "GGGATTAA");
        assert_eq!(str2, "GGTAA");
        println!("Alignment 1: {} alignment 2: {}", str1align, str2align);

        assert_eq!(str1align, "GGGATTAA");
        assert_eq!(str2align, "-GG--TAA");
    }


    #[test]
    fn test_large_alignment() {
        let scores = Scores::default_scores();

        let size = 50;
        let str1 = vec!['A'; size];
        let str2 = vec!['A'; size];

        for _ in 0..1000 {
            let _alignment = needleman_wunsch(&str1, &str2, &scores);
        }

        // assert_eq!(alignment.score, (size as f64)  * scores.match_score);
    }
}
