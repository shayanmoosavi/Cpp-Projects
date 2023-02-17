#pragma once
#include "matrix.h"


// Linear Algebra Namespace
namespace LinAlg {
	// Defining matrix to vector conversion
    /*Converts a column Matrix to a Vector.*/
    [[maybe_unused]] Vector ConvertToMatrix(const Matrix& A);

	// Creating an Identity Matrix with size dxd
    /*An Identity Matrix is defined as a Matrix which has zero elements everywhere except at the diagonal entries.*/
	Matrix Identity(int d);

	// Creating a Diagonal Matrix with diagonal elements
    /*A Diagonal Matrix is defined as a Matrix which has zero elements everywhere and at least one non-zero element
     at the diagonal entries.*/
    [[maybe_unused]] Matrix Diagonal(int d, const double var[]);

	// overloading the "*" operator for matrix scaler multiplication
	Matrix operator*(double c, const Matrix& A);

	// Adds the two matrices
    [[maybe_unused]] Matrix MatSum(const Matrix& A, const Matrix& B);

	// Calculates the matrix multiplication
    /*An "mxn" Matrix multiplied by an "nxp" Matrix is an "mxp" Matrix. Matrix Multiplication is
     undefined otherwise.*/
	Matrix MatMult(const Matrix& A, const Matrix& B);

	// Calculates the scaler c multiplied by matrix A
	Matrix MatScalerMult(double c, const Matrix& A);

	// Multiplies a row by a number
    [[maybe_unused]] Matrix RowMult(const Matrix& A, int row, double num);

	// Swaps two rows
    [[maybe_unused]] Matrix RowSwap(const Matrix& A, int row1, int row2);

	// Adds a multiple of row2 to row1
    [[maybe_unused]] Matrix RowAdd(const Matrix& A, int row1, int row2, double num);

	// Calculates the inverse of a matrix
    [[maybe_unused]] Matrix Inv(Matrix& A);

	// Transposing a vector
    /*Transposing a Vector yields a row matrix.*/
	Matrix VecTranspose(const Vector& v);

	// Overloading the "*" operator for vector scaler multiplication
	Vector operator*(double c, const Vector& v);

	// Adds the two vectors
	Vector Vecadd(const Vector& u, const Vector& v);

	// Subtracts the two vectors
	Vector Vecsub(const Vector& u, const Vector& v);

	// Multiplies a vector by a scaler
	Vector VecScalerMult(double c, const Vector& v);

	// Calculates the cross product
    [[maybe_unused]] Vector Cross(const Vector& u, const Vector& v);

	// Calculates the dot product
	double Dot(const Vector& u, const Vector& v);

	/*Solves the linear system Ax=b where A is a matrix,
	x is an unknown vector and b is a vector*/
    [[maybe_unused]] Matrix LinearSystemSolve(Matrix& A, const Vector& b);

	// Projects a vector v1 on the line spanned by v2
	Vector Proj(const Vector& v1, const Vector& v2);

	/*Applies the Gram-Schmidt process on the matrix A
	where columns of A are the set of vectors that needs to be made orthonormal.*/
    [[maybe_unused]] Matrix GramSchmidt(const Matrix& A);
}