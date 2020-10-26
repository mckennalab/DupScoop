mod convex;
mod kmer_orientation;
mod score_matrix;
pub mod mymatrix;
pub mod needleman;

extern crate bio;
extern crate clap;
extern crate csv;

extern crate matrix;

use std::str;
use std::fs::File;
use std::io::prelude::*;
use std::io::{Write};
use needleman::{Scores};
use clap::{Arg, App};

use bio::io::fastq;
use kmer_orientation::{ReferenceKmers, ReadOrientation};

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

    let kmer_orientation = ReferenceKmers::generate_kmers(&reference,&20);


    // read in sequences and align them
    // ----------------------------------------------------------------
    let scores = Scores::default_scores();
    let reader = fastq::Reader::from_file(input_file).unwrap();
    let mut output = File::create(output_file).unwrap();

    let mut count = 0;
    for result in reader.records() {
        // obtain record or fail with error
        let record = result.unwrap();

        str::from_utf8(record.seq()).unwrap().to_string();
        let orientation = kmer_orientation.vote_orientation(&str::from_utf8(record.seq()).unwrap(), &0.8, &50);

        let read_as_chars = match orientation {
            ReadOrientation::REV => {Some(ReferenceKmers::reverse_complement_sequence(str::from_utf8(record.seq()).unwrap()))}
            ReadOrientation::FWD => {Some(str::from_utf8(record.seq()).unwrap().to_string())}
            ReadOrientation::UNKNOWN => {None}
        };


        match read_as_chars {
            Some(p) => {
                let alignment = needleman::needleman_wunsch(&reference_as_chars, &p.chars().collect(), &scores);
                //convex::convex(&reference_as_chars, &p.chars().collect(), &scores);

                let str1align: String = alignment.seq_one_aligned.into_iter().collect();
                let str2align: String = alignment.seq_two_aligned.into_iter().collect();

                writeln!(output, ">{}", "ref")?;
                writeln!(output, "{}", str1align)?;
                writeln!(output, ">{}", record.id())?;
                writeln!(output, "{}", str2align)?;

                count += 1;
                if count % 1000 == 0 {
                    println!("Processed {} reads", count);
                }
            }
            None => {

            }
        }

    }
    //output.close();

    // we're ok!
    Ok(())
}

