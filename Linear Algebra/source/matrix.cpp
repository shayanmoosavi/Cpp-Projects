#include "matrix.h"
#include <iostream>
#include <iomanip>

namespace LinAlg {
	// EXCEPTIONS
	MatrixExceptions::MatrixExceptions(const char* msg) : message(msg) {}

	const char* MatrixExceptions::what() const throw() {
		return message;
	}

	// Matrix Constructor
	Matrix::Matrix(int row, int col) : Row(row), Col(col) {

		// Dynamically Allocating Memory
		Elem = new double* [row];
		for (int i = 0; i < row; i++) {
			Elem[i] = new double[col];
		}

	}

	// Matrix Destructor (Deleted for now)
	/*Matrix::~Matrix() {
		for (int i = 0; i < Row; i++) {
			delete Elem[i];
		}
		delete[] Elem;
		Elem = NULL;
	}*/

	// Defining vector to matrix conversion
	Matrix::Matrix(const Vector& v) {

		// Initializing the Rows and columns
		Row = v.Row;
		Col = v.Col;

		// Dynamically Allocating Memory
		Elem = new double* [Row];
		for (int i = 0; i < Row; i++) {
			Elem[i] = new double[Col];
		}
		// setting the elements
		for (int i = 0; i < Row; i++) {
			Elem[i][0] = v.Elem[i];
		}
	}

	// Asking the user to enter matrix elements
	void Matrix::EnterElements() {

		// Entering Matrix Elemnts
		for (int i = 0; i < Col; i++) {
			std::cout << "Enter col " << (i + 1) << ":\n";
			for (int j = 0; j < Row; j++) {
				std::cin >> Elem[j][i];
			}
		}
	}

	// Displaying the matrix
	void Matrix::Display() const {

		// Setting the precision for decimal places
		std::cout << std::setprecision(6) << std::fixed;

		// Displaying the matrix elements
		for (int i = 0; i < Row; i++) {
			for (int j = 0; j < Col; j++) {
				std::cout << Elem[i][j] << " ";
			}
			std::cout << std::endl;
		}
	}

	// Displaying the size of a matrix
	void Matrix::Size() const {
		std::cout << Row << 'x' << Col << std::endl;
	}

	// Overloading the "-" operator for matrix subtraction
	Matrix Matrix::operator-(const Matrix& A) {
		// Matrix dimensions should match
		if (Row != A.Row || Col != A.Col) {

			// Producing an error
			MatrixExceptions sub_err("Undefined Operation! Matrix dimensions do not match.");
			throw sub_err;
		}

		// Creating a new matrix object
		Matrix AminusB(A.Row, A.Col);

		// Doing the matrix addition
		for (int i = 0; i < A.Row; i++) {
			for (int j = 0; j < A.Col; j++) {
				AminusB.Elem[i][j] = Elem[i][j] - A.Elem[i][j];
			}
		}

		return AminusB;
	}

	// Calculating the trace of a matrix
	double Matrix::Trace() {

		// trace of a non-square matrix is undefined
		if (Row != Col) {

			// producing an error
			MatrixExceptions trace_err("Undefined operation! Trace is undefined for non-square matrices.");
			throw trace_err;
		}
		double sum = 0;
		for (int i = 0; i < Col; i++) {
			sum += Elem[i][i];
		}

		return sum;

	}

	// Calsulating the determinant of a matrix
	double Matrix::Det() {

		if (Row != Col) {
			MatrixExceptions det_err("ERROR: Determinant is undefined for non-square matrices.");
			throw det_err;
		}
		else {
			if (Row == 2) {
				return Elem[0][0] * Elem[1][1] - Elem[0][1] * Elem[1][0];
			}
			else {
				double sum = 0;
				for (int i = 0; i < Col; i++) {
					sum += pow(-1, i) * Elem[0][i] * SubMatrix(0, i).Det();
				}
				return sum;
			}
		}
	}

	// Transposing a matrix
	Matrix Matrix::Transpose() {

		// Creating a new matrix object with row and column numbers reversed
		Matrix TransMat(Col, Row);

		// Doing the transpose
		for (int i = 0; i < TransMat.Row; i++) {
			for (int j = 0; j < TransMat.Col; j++) {
				TransMat.Elem[i][j] = Elem[j][i];
			}
		}

		return TransMat;
	}

	// Creating a submatrix with row i and col j removed
	Matrix Matrix::SubMatrix(int row, int col) {
		int SubRow = Row - 1;
		int SubCol = Col - 1;
		Matrix sub(SubRow, SubCol);

		for (int i = 0; i < Row; i++) {
			if (i == row) {
				// Skipping the row
				continue;
			}
			else {
				for (int j = 0; j < Col; j++) {
					if (j == col) {
						// Skipping the col
						continue;
					}
					else {
						// Inserting the new matrix elements
						if (i > row) {
							if (j > col) {
								sub.Elem[i - 1][j - 1] = Elem[i][j];
							}
							else {
								sub.Elem[i - 1][j] = Elem[i][j];
							}
						}
						else {
							if (j > col) {
								sub.Elem[i][j - 1] = Elem[i][j];
							}
							else {
								sub.Elem[i][j] = Elem[i][j];
							}
						}
					}
				}
			}
		}

		return sub;
	}

	// Augment a matrix on the right with another matrix
	Matrix Matrix::MatAugment(const Matrix& A) {
		if (Row != A.Row) {
			MatrixExceptions aug_err("ERROR: Rows do not match.");
			throw aug_err;
		}
		else {
			Matrix AugMat(Row, Col + A.Col);
			for (int i = 0; i < AugMat.Row; i++) {
				for (int j = 0; j < AugMat.Col; j++) {
					if (j < Row) {
						AugMat.Elem[i][j] = Elem[i][j];
					}
					else {
						AugMat.Elem[i][j] = A.Elem[i][j - Col];
					}
				}
			}

			return AugMat;
		}

	}

	// Splitting a matrix into two other matrices sliced at specified column
	void Matrix::Split(int col, Matrix* Mat1, Matrix* Mat2) {
		if (col == 0) {
			throw MatrixExceptions("ERROR: Column index should be at least 1!");
		}
		else if (col == Col) {
			throw MatrixExceptions("ERROR: Cannot split at the last column!");
		}
		else {
			*Mat1 = Matrix(Row, col);
			*Mat2 = Matrix(Row, Col - col);
			for (int i = 0; i < Row; i++) {
				for (int j = 0; j < Col; j++) {
					if (j < col) {
						Mat1->Elem[i][j] = Elem[i][j];
					}
					else {
						Mat2->Elem[i][j - col] = Elem[i][j];
					}
				}
			}
		}
	}

	// Getting all the columns of a matrix
	void Matrix::GetCols() const {

		std::cout << std::setprecision(6) << std::fixed;
		for (int i = 0; i < Col; i++) {
			std::cout << "col " << i + 1 << std::endl;
			for (int j = 0; j < Row; j++) {
				std::cout << Elem[j][i] << std::endl;
			}
		}
	}

	// Checking whether the matrix is symmetric
	bool Matrix::isSymmetric() {

		// Symmetry is only defined for square matrices 
		if (Col != Row) {
			// producing an error
			MatrixExceptions sym_err("Symmetry is undefined for non-square matrices!");
			throw sym_err;
		}
		else {

			// checking the condition for symmetry (a_ij=a_ji)
			bool test = true;
			for (int i = 0; i < Row; i++) {
				for (int j = 0; j < Col; j++) {
					if (i == j) {
						continue;
					}
					else if (Elem[i][j] == Elem[j][i]) {
						test = true;
					}
					else {
						test = false;
					}
				}
			}
			if (test == true) {
				return true;
			}
			else {
				return false;
			}
		}
	}

	// Checking whether a matrix is a square matrix
	bool Matrix::isSquare() {
		if (Row == Col) {
			return true;
		}
		else {
			return false;
		}
	}
}