#pragma once
#include "vector.h"

namespace LinAlg {
	// EXCEPTIONS
	class MatrixExceptions : public std::exception {
		const char* message;
	public:

		MatrixExceptions(const char* msg);

		const char* what() const throw();
	};

	// Matrix Class
	class Matrix {
	public:
		int Row; // Row Number
		int Col; // Column Number
		double** Elem; // Matrix Elements

		// Matrix Constructor
		Matrix(int row, int col);

		// Matrix Destructor
		//~Matrix();

		// Defining vector to matrix conversion
		Matrix(const Vector& v);

		// Asking the user to enter matrix elements
		void EnterElements();

		// Displaying the matrix
		void Display() const;

		// Displaying the size of a matrix
		void Size() const;

		// Overloading the "+" operator for matrix addition
		Matrix operator+(const Matrix& A);

		// Overloading the "-" operator for matrix subtraction
		Matrix operator-(const Matrix& A);

		// Calculating the trace of a matrix
		double Trace();

		// Calsulating the determinant of a matrix
		double Det();

		// Transposing a matrix
		Matrix Transpose();

		// Creating a submatrix with row i and col j removed
		Matrix SubMatrix(int row, int col);

		// Augment a matrix on the right with another matrix
		Matrix MatAugment(const Matrix& A);

		// Splitting a matrix into two other matrices sliced at specified column
		void Split(int col, Matrix* Mat1, Matrix* Mat2);

		// Getting all the columns of a matrix
		void GetCols() const;

		// Checking whether the matrix is symmetric
		bool isSymmetric();

		// Checking whether a matrix is a square matrix
		bool isSquare();
	};
}