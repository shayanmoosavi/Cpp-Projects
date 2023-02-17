#pragma once
#include "vector.h"


// Linear Algebra Namespace
namespace LinAlg {
    // Custom exception class for vectors
    class MatrixExceptions : public std::exception {

        // Error Type
        //Represents what type of error caused an exception.
        const char *ErrorType;
        // Message
        // Explains the cause of the error.
        const char *Message;
    public:

        // The constructor of the MatrixExceptions class
        [[maybe_unused]] explicit MatrixExceptions(const char* errtyp, const char* msg);

        // Overriding the "what()" method of the base exceptions class
		[[nodiscard]] const char* what() const noexcept override;
	};

	// Matrix Class
	class Matrix {
	public:
		int Row; // Row Number
		int Col; // Column Number
		double** Elem; // Matrix Elements

		// The constructor for the Matrix class
        // Creates a Matrix object with the given row and column as arguments.
		Matrix(int row, int col);

		// The destructor for the Matrix class
        // Deletes the Matrix object.
		~Matrix();

		// Defining Vector to Matrix conversion
        // Converts a Vector to a Matrix
		explicit Matrix(const Vector& v);

		// Asks the user to enter Matrix elements
        [[maybe_unused]] void EnterElements() const;

		// Displays the Matrix
        [[maybe_unused]] void Display() const;

		// Displays the size of a Matrix
        [[maybe_unused]] void Size() const;

		// Overloading the "+" operator for Matrix Addition
		Matrix operator+(const Matrix& A) const;

		// Overloading the "-" operator for Matrix Subtraction
		Matrix operator-(const Matrix& A) const;

        // Defining the copy assignment operator
        // Copies the Matrix and assigns it to the left-hand side variable
        Matrix& operator=(const Matrix& A);

		// Calculates the trace of a Matrix.
        /*Trace of a Matrix is defined as the sum of the diagonal elements.
          Trace is only defined for square matrices.*/
        [[maybe_unused]] [[nodiscard]] double Trace() const;

		// Calculates the determinant of a Matrix.
        /*The determinant of a 2x2 Matrix of the form:
          a   b
          c   d*/
        /*is defined as ad-bc. The determinant of the higher dimension matrices are defined recursively.*/
		[[nodiscard]] double Det() const;

		// Transposing a Matrix
        /* Transpose of a Matrix is defined as exchanging the role of rows and columns. Meaning, rows of the
         Matrix are the columns of the transposed Matrix and vise-versa.*/
        [[maybe_unused]] [[nodiscard]] Matrix Transpose() const;

		// Creating a sub-matrix with row i and column j removed
        /* A sub-matrix is defined as the Matrix that we get if we ignore row i and column j.*/
		[[nodiscard]] Matrix SubMatrix(int row, int col) const;

		// Augments a Matrix on the right with another matrix
        /* Augmenting a Matrix is defined as appending two matrices together.*/
		[[nodiscard]] Matrix MatAugment(const Matrix& A) const;

		// Splitting a matrix into two other matrices sliced at specified column
		void Split(int col, Matrix* Mat1, Matrix* Mat2) const;

		// Getting all the columns of a matrix
        /* Displays all the columns of a Matrix in the terminal.*/
        [[maybe_unused]] void GetCols() const;

		// Checking whether the matrix is symmetric
        /* Returns "true" if the Matrix is symmetric and "false" otherwise.*/
        [[maybe_unused]] [[nodiscard]] bool isSymmetric() const;

		// Checking whether a matrix is a square matrix
        /* Returns "true" if the Matrix is a square Matrix and "false" otherwise.*/
        [[maybe_unused]] [[nodiscard]] bool isSquare() const;
	};
}