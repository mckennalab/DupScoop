use std::fmt::{Display, Debug};
use std::borrow::Cow;

/// A row-first matrix class
pub struct Matrix<T: Clone> {
    values: Vec<T>,
    row_length: usize,
}

impl<T> Matrix<T> where T: Clone + Debug + Sized {
    #[inline]
    pub fn cols(&self) -> usize {
        self.row_length
    }

    pub fn rows(&self) -> usize {
        self.values.len() / self.row_length
    }
    
    #[inline]
    pub fn get<'a>(&self, row: usize, col: usize) -> T {
        self.values[(col * self.row_length) + row].clone()
    }
    
    #[inline]
    pub fn set(&mut self, row: usize, col: usize, value: T) {
        self.values[(col * self.row_length) + row] = value
    }

    pub fn new(row: usize, col: usize, initialize: T) -> Matrix<T> {
        let mat = Matrix::<T>{
            values: vec![initialize; row * col],
            row_length: row,
        };
        mat
    }

    pub fn row(&self, row: usize) -> Vec<T> {
        self.values[(row * self.row_length)..((row + 1) * self.row_length)].to_vec()
    }
    
    pub fn col(&self,col: usize) -> Vec<T> {
        let mut mt = Vec::new();
        for i in 0..self.rows() {
            mt.push(self.get(i,col))
        };
        mt
    }

    pub fn print_matrix(&self) {
        for ix in 0..self.rows() {
            for iy in 0..self.cols() {
                print!("{:?}",self.get(ix,iy));
                print!("{}",",");
            }
            print!("{}","\n");
        }
    }
}
