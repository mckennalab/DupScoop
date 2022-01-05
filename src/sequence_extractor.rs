use std::f64;
use std::fmt;
use string_builder::Builder;
/*
// an extracted sequence from a read/reference alignment
pub struct ReferenceExtract {
    reference: String,
    read: String,
    ref_start: usize,
    length: usize,
}

impl ReferenceExtract {

    #[inline]
    pub fn length(&self) -> usize {self.length.clone()}

    #[inline]
    pub fn ref_start(&self) -> usize {self.ref_start.clone()}

    #[inline]
    pub fn read(&self) -> String {self.read.clone()}

    #[inline]
    pub fn reference(&self) -> String {self.reference.clone()}
}


pub fn extract_alignment(reference: String, read: String, start: usize, length: usize) -> ReferenceExtract {
    let mut ref_position : usize = 0;
    let mut offset : usize = 0;

    let mut stop = start + length;

    let mut ref_start = usize::MAX;
    let mut ref_stop = 0;

    let ref_as_vec = reference[offset].chars().to_vec();

    while ref_position < stop {
        match ref_as_vec(ref_position).unwrap() {
            _x if offset < start => {},
            '-' => {
                offset += 1;
            }
            _x if ref_position < start => {

            },
            '-' => {
                offset += 1;
            }
            _ => {
                reference_builder.append(reference[offset].to_vec().unwrap());
                read_builder.append(read[offset].to_vec().unwrap());
                offset += 1;
                ref_position += 1;
            }
        }
    }

    ReferenceExtract {
        reference: reference_builder.string().unwrap(),
        read: read_builder.string().unwrap(),
        ref_start: start.clone(),
        length: length.clone(),
    }
}
*/