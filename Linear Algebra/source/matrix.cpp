#include "matrix.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstring>


// Linear Algebra Namespace
namespace LinAlg {

    [[maybe_unused]] MatrixExceptions::MatrixExceptions(
            // Member Initializer List
            const char* errtyp, const char* msg) : ErrorType(errtyp), Message(msg) {}

	const char* MatrixExceptions::what() const noexcept {
        // Concatenating the strings
        static char buffer[sizeof(this->ErrorType) + sizeof(" - ") + sizeof(this->Message)];
        strcpy(buffer, ErrorType);
        strcat(buffer, " - ");
        strcat(buffer, Message);

        return buffer;
	}

	Matrix::Matrix(int row, int col) : Row(row), Col(col) {

		// Dynamically Allocating Memory
		Elem = new double* [row];
		for (int i = 0; i < row; i++) {
			Elem[i] = new double[col];
		}

	}

	Matrix::~Matrix() {
		for (int i = 0; i < Row; i++) {
			delete Elem[i];
		}
		delete[] Elem;
	}

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

    [[maybe_unused]] void Matrix::EnterElements() const {

		// Entering Matrix Elements
		for (int i = 0; i < Col; i++) {
			std::cout << "Enter col " << (i + 1) << ":\n";
			for (int j = 0; j < Row; j++) {
				std::cin >> Elem[j][i];
			}
		}
	}

    [[maybe_unused]] void Matrix::Display() const {

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

    [[maybe_unused]] void Matrix::Size() const {
		std::cout << Row << 'x' << Col << std::endl;
	}

    Matrix Matrix::operator+(const Matrix& A) const {

        // Matrix dimensions should match
        if (Row != A.Row || Col != A.Col) {

            // Producing an error if the dimensions do not match
            throw MatrixExceptions("[NonMatchingDimensionsException]",
                                   "Undefined Operation! Matrix dimensions do not match.");

        }

        /* Creating a new matrix object. Declared static in order to prevent the destructor
         from being called inappropriately.*/
        static Matrix AplusB(A.Row, A.Col);

        // Doing the matrix addition
        for (int i = 0; i < A.Row; i++) {
            for (int j = 0; j < A.Col; j++) {
                AplusB.Elem[i][j] = Elem[i][j] + A.Elem[i][j];
            }
        }

        return AplusB;

    }


    Matrix Matrix::operator-(const Matrix& A) const {

        // Matrix dimensions should match
		if (Row != A.Row || Col != A.Col) {
			// Producing an error message if Matrix dimensions do not match
			throw MatrixExceptions("[NonMatchingDimensionsException]",
                                   "Undefined Operation! Matrix dimensions do not match.");
		}

		/* Creating a new matrix object. Declared static in order to prevent the destructor
         from being called inappropriately.*/
		static Matrix AminusB(A.Row, A.Col);

		// Doing the matrix subtractions
		for (int i = 0; i < A.Row; i++) {
			for (int j = 0; j < A.Col; j++) {
				AminusB.Elem[i][j] = Elem[i][j] - A.Elem[i][j];
			}
		}

		return AminusB;
	}

    Matrix &Matrix::operator=(const Matrix &A) {

        if(this != &A){ // checking for self assignment

            this->Row = A.Row;
            this->Col = A.Col;

            for (int i = 0; i < this->Row; i++) {
                // Assigning the values to the left-hand side variables
                for (int j = 0; j < this->Col; j++){
                    this->Elem[i][j] = A.Elem[i][j];
                }
            }
        }

        return *this;
    }

    [[maybe_unused]] double Matrix::Trace() const {

		// trace of a non-square matrix is undefined
		if (Row != Col) {
			// Producing an error message if it is a non-square Matrix
			throw MatrixExceptions("[UndefinedOperationException]",
                                   "Undefined operation! Trace is undefined for non-square matrices.");
		}

        // Summing the diagonal elements of the Matrix
		double sum = 0;
		for (int i = 0; i < Col; i++) {
			sum += Elem[i][i];
		}

		return sum;

	}

	double Matrix::Det() const {

        // The determinant of a non-square matrix is undefined
		if (Row != Col) {
            // Producing an error message if it is a non-square Matrix
            throw MatrixExceptions("[UndefinedOperationException]",
                                   "Undefined operation! Determinant is undefined for non-square matrices.");
		}
		else {

			if (Row == 2) {
                // Determinant of a 2x2 Matrix
				return Elem[0][0] * Elem[1][1] - Elem[0][1] * Elem[1][0];
			}
			else {
				// Recursive definition of the determinant
                double sum = 0;
				for (int i = 0; i < Col; i++) {
					sum += pow(-1, i) * Elem[0][i] * SubMatrix(0, i).Det();
				}

				return sum;
			}
		}
	}

    [[maybe_unused]] Matrix Matrix::Transpose() const {

		/* Creating a new matrix object with row and column numbers reversed. Declared static in order to
		 prevent the destructor from being called inappropriately*/
		static Matrix TransMat(Col, Row);

		// Doing the transpose
		for (int i = 0; i < TransMat.Row; i++) {
			for (int j = 0; j < TransMat.Col; j++) {
				TransMat.Elem[i][j] = Elem[j][i];
			}
		}

		return TransMat;
	}

	Matrix Matrix::SubMatrix(int row, int col) const {

        int SubRow = Row - 1;
		int SubCol = Col - 1;

        static Matrix sub(SubRow, SubCol);

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

	Matrix Matrix::MatAugment(const Matrix& A) const {

		if (Row != A.Row) {
            // Rows must match when we want to append a Matrix to the right
			throw MatrixExceptions("[NonMatchingDimensionsException]","ERROR: Rows do not match.");
		}
		else {

			static Matrix AugMat(Row, Col + A.Col);

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

	void Matrix::Split(int col, Matrix* Mat1, Matrix* Mat2) const {

		if (col == 0) {
            // At least one column is needed for splitting into two matrices
            throw MatrixExceptions("[MiscException]","ERROR: Column index should be at least 1!");
		}
		else if (col == Col) {
            // At least one column is needed for splitting into two matrices
            throw MatrixExceptions("[MiscException]","ERROR: Cannot split at the last column!");
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

    [[maybe_unused]] void Matrix::GetCols() const {

		std::cout << std::setprecision(6) << std::fixed;
		for (int i = 0; i < Col; i++) {
			std::cout << "col " << i + 1 << std::endl;
			for (int j = 0; j < Row; j++) {
				std::cout << Elem[j][i] << std::endl;
			}
		}
	}

    [[maybe_unused]] bool Matrix::isSymmetric() const {

        // Non-square matrices are not symmetric
		if (Col != Row) {
			return false;
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
			if (test) {
				return true;
			}
			else {
				return false;
			}
		}
	}

    [[maybe_unused]] bool Matrix::isSquare() const {
		if (Row == Col) {
			return true;
		}
		else {
			return false;
		}
	}
}