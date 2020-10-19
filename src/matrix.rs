use std::fmt::{Display, Debug};
use std::borrow::Cow;
use std::ops;

use std::f64;

/// A row-first matrix class
pub struct Matrix<T: Clone> {
    values: Vec<T>,
    row_length: usize,
}

impl<T> Matrix<T> where T: Clone + Debug + Sized {
    #[inline]
    pub fn rows(&self) -> usize {
        self.row_length
    }

    pub fn cols(&self) -> usize {
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

    #[inline]
    pub fn set_column_range(&mut self, row_start: usize, row_end: usize, col: usize, value: T) {
        for row_index in row_start..row_end {
            self.values[(col * self.row_length) + row_index] = value.clone();
        }
    }

    #[inline]
    pub fn set_row_range(&mut self, row: usize, col_start: usize, col_end: usize, value: T) {
        for col_index in col_start..col_end {
            self.values[(col_index * self.row_length) + row] = value.clone();
        }
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


pub fn maximize_over_column(mtx: &Matrix<f64>, col: usize, cur_row: usize, scoring_function: &Fn(usize, f64) -> f64) -> (usize, f64) {

    let mut max_index = cur_row;
    let mut max_score = f64::MIN;
    for n in (0..cur_row).rev() {
        let score = scoring_function(cur_row - n, mtx.get(n,col));
        if score > max_score {
            max_index = n;
            max_score = score;
        }

    };
    (max_index,max_score)
}

pub fn maximize_over_row(mtx: &Matrix<f64>, row: usize, cur_column: usize, scoring_function: &Fn(usize, f64) -> f64) -> (usize, f64) {
    let mut max_index = cur_column;
    let mut max_score = f64::MIN;
    for n in (0..cur_column).rev() {
        let score = scoring_function(cur_column - n, mtx.get(row,n));
        if score > max_score {
            max_index = n;
            max_score = score;
        }

    };
    (max_index,max_score)
}


#[cfg(test)]
mod tests {
    // Note this useful idiom: importing names from outer (for mod tests) scope.
    use super::*;

    #[test]
    fn basic_size_setup() {

        let mut mtx = Matrix::new(10, 5, 0.0);

        assert_eq!(mtx.rows(), 10);
        assert_eq!(mtx.cols(), 5);
    }

    #[test]
    fn get_set() {

        let mut mtx = Matrix::new(10, 5, 0.0);

        mtx.set(1,4, 2.0);
        assert_eq!(mtx.get(1,4),2.0);
        assert_eq!(mtx.get(0,4),0.0);
    }

    #[test]
    fn minimize_over_x_test() {

        let mut mtx = Matrix::new(10, 5, 0.0);
        let convert_function = |distance: usize, score_at_distance: f64| (-10.0 - ((distance as f64) * 0.1)) + score_at_distance;

        mtx.set(3,4, 20.0);
        let (index, score) = maximize_over_column(&mtx, 4,4, &convert_function);
        assert_eq!(score,9.9);
        let (index, score) = maximize_over_row(&mtx, 4,4, &convert_function);
        assert_eq!(score,-10.1);

    }
}
