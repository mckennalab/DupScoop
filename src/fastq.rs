use std::fs::File;
use std::io::{Write, BufReader, BufRead};

#[derive(Debug)]
pub struct Fastq {
    read_name: Vec<u8>,
    bases: Vec<u8>,
    quals: Vec<u8>,
}
#[derive(Debug)]
pub struct FastqString {
    read_name: Vec<u8>,
    bases: Vec<u8>,
    quals: Vec<u8>,
}


pub fn file_input(path: String) -> Fastq {
    let f = File::open(path).expect("Unable to open file");
    let reader = BufReader::new(f);
    
    let mut line1 = String::new();
    
    try!(reader.read_line(&mut line1));

    let mut line2 = String::new();
    try!(reader.read_line(&mut line2));

    let mut line3 = String::new();
    try!(reader.read_line(&mut line3));

    let mut line4 = String::new();
    try!(reader.read_line(&mut line4));
    
    return Fastq{
        read_name: line1.as_bytes().to_vec(),
        bases: line2.as_bytes().to_vec(),
        quals: line4.as_bytes().to_vec(),
    }
}


#[cfg(test)]
mod tests {
    // Note this useful idiom: importing names from outer (for mod tests) scope.
    use super::*;

    #[test]
    fn test_basic_alignment() {
        let filepath = "/Users/aaronmck/Desktop/code/SCIMaul/data/cells/cell.index1.fq";
        file_input(&filepath);
        //assert_eq!(alignment.score, 3.0 * scores.match_score);
    }
}
