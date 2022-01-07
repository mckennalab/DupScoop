use mymatrix;
use needleman::Direction::{Left, Up, Diag, Done};
use std::f64;
use std::fmt;

#[allow(dead_code)]
pub struct Scores {
    pub match_score: f64,
    pub mismatch_score: f64,
    pub gap_open: f64,
    pub gap_ext: f64,
    pub gap_start: f64,
    pub gap_end: f64,
}

impl Scores {
    pub fn default_scores() -> Scores {
        Scores {
            match_score: 6.0,
            mismatch_score: -5.0,
            gap_open: -10.0,
            gap_ext: -6.0,
            gap_start: -10.0,
            gap_end: -10.0,
        }
    }

    pub fn plasmid_aligmment_scores() -> Scores {
        Scores {
            match_score: 2.0,
            mismatch_score: -6.0,
            gap_open: -10.0,
            gap_ext: -0.5,
            gap_start: -10.0,
            gap_end: -10.0,
        }
    }

    #[inline]
    pub fn scoring_function(base1: char, base2: char, scores: &Scores) -> f64 {
        if base1 == base2 {
            return scores.match_score;
        } else if base1 == 'N' || base2 == 'N' || base1 == 'Y' || base2 == 'Y' {
            return 0.0;
        } else {
            return scores.mismatch_score;
        }
    }
}


#[allow(dead_code)]
#[derive(Copy, Clone, Debug, PartialEq)]
pub enum Direction {
    Up,
    Left,
    Diag,
    Done,
}

impl fmt::Display for Direction {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{:?}", self)
    }
}

#[allow(dead_code)]
pub struct Alignment {
    pub seq_one: Vec<char>,
    pub seq_two: Vec<char>,
    pub score: f64,
    pub start_x: usize,
    pub start_y: usize,
    pub end_x: usize,
    pub end_y: usize,
    pub seq_one_aligned: Vec<char>,
    pub seq_two_aligned: Vec<char>,
}

/// Aligns two sequences using the Needleman Wunsch global alignment with simple gap scoring 
pub fn needleman_wunsch(seq1: &Vec<char>, seq2: &Vec<char>, scores: &Scores) -> Alignment {
    let seq1_limit = seq1.len() + 1;
    let seq2_limit = seq2.len() + 1;

    let mut mtx = mymatrix::MyMatrix::new(seq1_limit, seq2_limit, 0.0);
    println!("made matrix of {} {}", mtx.rows(), mtx.cols());
    let mut trc: mymatrix::MyMatrix<Direction> = mymatrix::MyMatrix::new(seq1_limit, seq2_limit, Done);
    needleman_wunsch_borrow(seq1, seq2, &mut mtx, &mut trc, scores)
}

pub fn needleman_wunsch_borrow(seq1: &Vec<char>,
                               seq2: &Vec<char>,
                               mtx: &mut mymatrix::MyMatrix<f64>,
                               trc: &mut mymatrix::MyMatrix<Direction>,
                               scores: &Scores) -> Alignment {

    let seq1_limit = seq1.len() + 1;
    let seq2_limit = seq2.len() + 1;
    // first square
    mtx.set(0, 0, 0.0);

    // initialize the top row and first column
    for n in 1..seq1_limit {
        mtx.set(n, 0, scores.gap_ext * (n as f64));
        trc.set(n, 0, Up);
    }
    for n in 1..seq2_limit {
        mtx.set(0, n, scores.gap_ext * (n as f64));
        trc.set(0, n, Left);
    }


    // fill in the matrix
    for ix in 1..seq1_limit {
        for iy in 1..seq2_limit {
            let score = Scores::scoring_function(seq1[ix - 1], seq2[iy - 1], scores);
            if (ix == iy || ((iy % seq1.len()) == ix) || ((ix % seq2.len()) == iy)) && ix % 100 == 0 {
                println!("NOGO is from {},{}", ix, iy);
            }
            let up_t = (mtx.get(ix - 1, iy) + scores.gap_ext, Up);
            let left_t = (mtx.get(ix, iy - 1) + scores.gap_ext, Left);
            let diag_t = (mtx.get(ix - 1, iy - 1) + score, Diag);

            //let max = max2(max2(max2(upT, leftT), diagT), (0.0,Diag));
            let max = max2(max2(up_t, left_t), diag_t);


            if ix % 1000 == 0 && iy % 1000 == 0 {
                println!("top is now {} from {},{}", max.0, ix, iy);
            }

            mtx.set(ix, iy, max.0);
            trc.set(ix, iy, max.1);

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

        let up = (30.0,Up);
        let left = (40.0,Left);
        let diag = (50.0,Diag);
        let max = max2(max2(diag, left), up);
        assert_eq!(max.0, 50.0);
        assert_eq!(max.1, Diag);

        let up = (60.0,Up);
        let left = (40.0,Left);
        let diag = (50.0,Diag);
        let max = max2(max2(diag, left), up);
        assert_eq!(max.0, 60.0);
        assert_eq!(max.1, Up);

        let up = (30.0,Up);
        let left = (60.0,Left);
        let diag = (50.0,Diag);
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
