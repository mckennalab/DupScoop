use matrix;

pub struct Scores {
    match_score: f64,
    mismatch_score: f64,
    gap_open: f64,
    gap_ext: f64,
    gap_start: f64,
    gap_end: f64,
}

impl Scores {
    pub fn default_scores() -> Scores {
        Scores {
            match_score: 4.0,
            mismatch_score: -5.0,
            gap_open: -10.0,
            gap_ext: -10.0,
            gap_start: -10.0,
            gap_end: -10.0,
        }
    }
}

#[derive(Copy, Clone, Debug)]
pub enum TracebackDirection {
    UP,
    LEFT,
    DIAG,
}


pub struct Alignment {
    seq_one: Vec<char>,
    seq_two: Vec<char>,
    score: f64,
    seq_one_aligned: Vec<char>,
    seq_two_aligned: Vec<char>,
}

pub struct AligmentStructure {
    mtx: matrix::Matrix<f64>,
    trc: matrix::Matrix<f64>,
}

/// Aligns two sequences using the Needleman Wunsch global alignment with simple gap scoring 
pub fn needleman_wunsch(seq1: &Vec<char>, seq2: &Vec<char>, scores: &Scores) -> Alignment {
    let seq1_limit = seq1.len() + 1;
    let seq2_limit = seq2.len() + 1;
    
    let mut mtx = matrix::Matrix::new(seq1_limit, seq2_limit, 0.0);
    let mut trc = matrix::Matrix::new(seq1_limit, seq2_limit, TracebackDirection::DIAG);
    
    // first square
    mtx.set(0,0,0.0);
    
    // initialize the top and first column
    for n in 1..seq1_limit {
        mtx.set(0,n,scores.gap_open + (scores.gap_ext * (n as f64)));
        mtx.set(n,0,scores.gap_open + (scores.gap_ext * (n as f64)));
    }
    mtx.print_matrix();

    // fill in the matrix
    for ix in 1..seq1_limit {
        for iy in 1..seq2_limit {
            let score = if seq1[ix - 1] == seq2[iy - 1] { scores.match_score } else { scores.mismatch_score };
            let up    = mtx.get(ix - 1,iy - 1) + scores.gap_ext;
            let left  = mtx.get(ix,iy - 1) + scores.gap_ext;
            let diag  = mtx.get(ix - 1,iy - 1) + score;
            
            // Find the minimum of d11, d01, d10
            // by enumerating all the cases. 
            if up < left {
                if diag < left {
                    mtx.set(ix,iy,left);
                    trc.set(ix,iy,TracebackDirection::LEFT);
                } else {
                    mtx.set(ix,iy,diag);
                    trc.set(ix,iy,TracebackDirection::DIAG);
                }
            } else {
                if diag < up {
                    mtx.set(ix,iy,up);
                    trc.set(ix,iy,TracebackDirection::UP);
                } else {
                    mtx.set(ix,iy,diag);
                    trc.set(ix,iy,TracebackDirection::DIAG);
                }
            };
        }
    }
    traceback(seq1,seq2,trc,mtx.get(seq1_limit-1,seq2_limit-1))
}

/// Aligns two sequences using the Convex alignment with an arbritrary scoring function (that should be convex of distance)
pub fn convex_alignment(seq1: &Vec<char>, seq2: &Vec<char>, scores: &Scores, scoring_function: &Fn(usize, f64) -> f64) -> Alignment {
    let seq1_limit = seq1.len() + 1;
    let seq2_limit = seq2.len() + 1;
    
    let mut mtx = matrix::Matrix::new(seq1_limit, seq2_limit, 0.0);
    let mut trc = matrix::Matrix::new(seq1_limit, seq2_limit, TracebackDirection::DIAG);
    
    // first square
    mtx.set(0,0,0.0);
    
    // initialize the top and first column
    for n in 1..seq1_limit {
        mtx.set(0,n,scores.gap_open + (scores.gap_ext * (n as f64)));
        mtx.set(n,0,scores.gap_open + (scores.gap_ext * (n as f64)));
    }
    mtx.print_matrix();

    // fill in the matrix
    for ix in 1..seq1_limit {
        for iy in 1..seq2_limit {
            let score = if seq1[ix - 1] == seq2[iy - 1] { scores.match_score } else { scores.mismatch_score };

            let row = mtx.row(ix);
            let col = mtx.col(iy);
                
            let up    = mtx.get(ix - 1,iy - 1) + scores.gap_ext;
            let left  = mtx.get(ix,iy - 1) + scores.gap_ext;
            let diag  = mtx.get(ix - 1,iy - 1) + score;
            
            // Find the minimum of d11, d01, d10
            // by enumerating all the cases. 
            if up < left {
                if diag < left {
                    mtx.set(ix,iy,left);
                    trc.set(ix,iy,TracebackDirection::LEFT);
                } else {
                    mtx.set(ix,iy,diag);
                    trc.set(ix,iy,TracebackDirection::DIAG);
                }
            } else {
                if diag < up {
                    mtx.set(ix,iy,up);
                    trc.set(ix,iy,TracebackDirection::UP);
                } else {
                    mtx.set(ix,iy,diag);
                    trc.set(ix,iy,TracebackDirection::DIAG);
                }
            };
        }
    }
    traceback(seq1,seq2,trc,mtx.get(seq1_limit-1,seq2_limit-1))
}


/// Produces an alignment object from a traceback matrix and two sequences
pub fn traceback(seq1: &Vec<char>, seq2: &Vec<char>, trc: matrix::Matrix<TracebackDirection>, score: f64) -> Alignment {
    assert!(seq1.len() == trc.rows() - 1, "The matrix doesn't have the right number of rows");
    assert!(seq2.len() == trc.cols() - 1, "The matrix doesn't have the right number of columns");

    let mut alignment1 = Vec::new();
    let mut alignment2 = Vec::new();

    let mut index1 = trc.rows() - 1;
    let mut index2 = trc.cols() - 1;

    let gap = '-';
    
    while index1 > 0 || index2 > 0 {
        match trc.get(index1,index2) {
            TracebackDirection::LEFT => {
                alignment1.push(seq1[index1 - 1]);
                alignment2.push(gap);
                index1 -= 1;
            }
            TracebackDirection::DIAG => {
                alignment1.push(seq1[index1 - 1]);
                alignment2.push(seq2[index2 - 1]);
                index1 -= 1;
                index2 -= 1;
            }
            TracebackDirection::UP => {
                alignment1.push(gap);
                alignment2.push(seq2[index2 - 1]);
                index2 -= 1;
            }
        }
    }
    
    while index1 > 0 {
        alignment1.push(seq1[index1 -1]);
        index1 -= 1;
    }
    while index2 > 0 {
        alignment2.push(seq2[index2 -1]);
        index2 -= 1;
    }

    alignment1.reverse();
    alignment2.reverse();
    
    return Alignment {
        seq_one: seq1.to_vec(),
        seq_two: seq2.to_vec(),
        score: score,
        seq_one_aligned: alignment1,
        seq_two_aligned: alignment2,
    }
}


#[cfg(test)]
mod tests {
    // Note this useful idiom: importing names from outer (for mod tests) scope.
    use super::*;

    #[test]
    fn test_basic_alignment() {
        let scores = Scores::default_scores();
        let alignment = needleman_wunsch(&vec!['A','A','A'], &vec!['A','A','A'], &scores);

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
    fn test_mismatch_alignment() {
        let scores = Scores::default_scores();
        let alignment = needleman_wunsch(&vec!['A','T','A'], &vec!['A','A','A'], &scores);

        let str1: String = alignment.seq_one.into_iter().collect();
        let str2: String = alignment.seq_two.into_iter().collect();
        let str1align: String = alignment.seq_one_aligned.into_iter().collect();
        let str2align: String = alignment.seq_two_aligned.into_iter().collect();
        
        assert_eq!(str1, "ATA");
        assert_eq!(str2, "AAA");
        assert_eq!(str1align, "ATA");
        assert_eq!(str2align, "AAA");
        assert_eq!(alignment.score, 2.0 * scores.match_score + scores.mismatch_score);
    }

    #[test]
    fn test_large_alignment() {
        let scores = Scores::default_scores();

        let size = 200;
        let str1 = vec!['A'; size];
        let str2 = vec!['A'; size];

        for i in 0..10000 {
            let alignment = needleman_wunsch(&str1, &str2, &scores);
        }

        // assert_eq!(alignment.score, (size as f64)  * scores.match_score);
    }
 }
