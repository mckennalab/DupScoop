mod kmer_orientation;
mod score_matrix;
pub mod matrix;
pub mod needleman;

extern crate bio;
extern crate clap;
extern crate csv;

use std::str;
use std::error::Error;
use std::io;
use std::process;
use std::fs::File;
use std::io::prelude::*;
use std::io::{Write, BufReader, BufRead};
use needleman::{Scores, Direction};
use clap::{Arg, App, SubCommand};

use bio::alphabets;
use bio::data_structures::suffix_array::suffix_array;
use bio::data_structures::bwt::{bwt, less, Occ};
use bio::data_structures::fmindex::{FMIndex, FMIndexable};
use bio::io::fastq;
use needleman::Direction::Done;

fn main() -> std::io::Result<()> {
    let matches = App::new("RustyShovel")
        .version("1.0")
        .author("Aaron M. <aaronmck@uw.edu>")
        .about("Various aligner implementations for digging around in genomic data")
        .arg(Arg::with_name("input")
             .short("i")
             .long("input")
             .value_name("FILE")
             .help("An input file containing input reads")
             .takes_value(true))
        .arg(Arg::with_name("reference")
            .short("r")
            .long("ref")
            .value_name("FILE")
            .help("The reference we will align to")
            .takes_value(true))
        .arg(Arg::with_name("output")
             .short("o")
             .long("output")
             .value_name("FILE")
             .help("the alignment output file")
             .takes_value(true))
        .arg(Arg::with_name("aligner")
             .short("a")
             .long("aligner")
             .value_name("FILE")
             .help("the alignment method")
             .takes_value(true))
        .arg(Arg::with_name("scoreMatrix")
            .short("s")
            .long("scores")
            .value_name("FILE")
            .help("a set of pairs detailing the cost/reward of matching base x to base y")
            .takes_value(true))
        .get_matches();

    // pull out the command line arguments and setup input and output
    // ----------------------------------------------------------------
    let input_file = matches.value_of("input").unwrap_or("input.fq");
    let output_file = matches.value_of("output").unwrap_or("output.fa");
    let reference_file = matches.value_of("reference").unwrap_or("reference.fa");
    let aligner = matches.value_of("aligner").unwrap_or("swa");
    println!("input : {}, output : {} aligner : {}",input_file, output_file, aligner);

    // setup the output fasta file
    let output = File::create(output_file)?;
    //write!(output, "Rust\n\nFun")?;

    // read in the reference, and covert to a string
    // ----------------------------------------------------------------
    let mut file = File::open(reference_file)?;
    let mut reference = String::new();
    file.read_to_string(&mut reference)?;

    // slice off the first line (name), and combine the rest into a single string of characters
    let first_line_marker = reference.find("\n");
    println!("first endline: {}" ,first_line_marker.unwrap());

    match first_line_marker {
        // The division was valid
        Some(x) => {
            reference = reference.split_off(x).replace("\n", "");
        }
        // The division was invalid
        None    => panic!("We couldn't find the reference name in your input file {}", input_file),
    }

    let reference_as_chars: Vec<char> = reference.to_string().chars().collect();

    // read in sequences and align them
    // ----------------------------------------------------------------
    let scores = Scores::default_scores();
    
    let reader = fastq::Reader::from_file(input_file).unwrap();

    let seq1_limit = reference_as_chars.len() + 1;

    let mut output = File::create(output_file).unwrap();

    let mut count = 0;
    for result in reader.records() {
        // obtain record or fail with error
        let record = result.unwrap();

        let read_as_chars: Vec<char> = str::from_utf8(record.seq()).unwrap().to_string().chars().collect();
        let seq2_limit = read_as_chars.len() + 1;

        //let mut mtx = matrix::Matrix::new(seq1_limit, seq2_limit, 0.0);
        //let mut trc: matrix::Matrix<Direction> = matrix::Matrix::new(seq1_limit, seq2_limit, Done);

        //let alignment = needleman::needleman_wunsch_borrow(&reference_as_chars, &read_as_chars,&mut mtx, &mut trc, &scores);
        let alignment = needleman::needleman_wunsch(&reference_as_chars, &read_as_chars, &scores);
        let str1: String = alignment.seq_one.into_iter().collect();
        let str2: String = alignment.seq_two.into_iter().collect();
        let str1align: String = alignment.seq_one_aligned.into_iter().collect();
        let str2align: String = alignment.seq_two_aligned.into_iter().collect();

        writeln!(output, ">{}", "ref")?;
        writeln!(output, "{}", str1align)?;
        writeln!(output, ">{}", record.id())?;
        writeln!(output, "{}", str2align)?;

        count += 1;
        if count % 50 == 0 {
            println!("Processed {} reads", count);
        }
    }
    //output.close();

    // we're ok!
    Ok(())
}

