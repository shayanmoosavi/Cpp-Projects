#pragma once
#include "matrix.h"

namespace LinAlg {
	// Defining matrix to vector conversion
	Vector ConvertToMatrix(const Matrix& A);

	// Creating an identity matrix with size dxd
	Matrix Identity(int d);

	// Creating a diagonal matrix with diagonal elements
	Matrix Diagonal(const int d, double var[]);

	// overloading the "*" operator for matrix scaler multiplication
	Matrix operator*(double c, const Matrix& A);

	// Calculating the matrix addition
	Matrix MatSum(const Matrix& A, const Matrix& B);

	// Calculating the matrix multiplication
	Matrix MatMult(const Matrix& A, const Matrix& B);

	// Calculating the scaler c multiplied by natrix A
	Matrix MatScalerMult(double c, const Matrix& A);

	// Multiplying a row by a number
	Matrix RowMult(const Matrix& A, int row, double num);

	// Swapping two rows
	Matrix RowSwap(const Matrix& A, int row1, int row2);

	// Adding row2 multiplied by a number to row1
	Matrix RowAdd(const Matrix& A, int row1, int row2, double num);

	// Calculating the inverse of a matrix
	Matrix Inv(Matrix& A);

	// Transposing a vector
	Matrix VecTranspose(const Vector& v);

	// Overloading the "*" operator for vector scaler multiplication
	Vector operator*(double c, const Vector& v);

	// Adding two vectors
	Vector Vecadd(const Vector& u, const Vector& v);

	// Subtracing two vectors
	Vector Vecsub(const Vector& u, const Vector& v);

	// Multiplying a vector by a scaler
	Vector VecScalerMult(double c, const Vector& v);

	// Calculating the cross product
	Vector Cross(const Vector& u, const Vector& v);

	// Calculating the dot product
	double Dot(const Vector& u, const Vector& v);

	/*Solves the linear system Ax=b where A is a matrix,
	x is an unknown vector and b is a vector*/
	Matrix LinearSystemSolve(Matrix& A, const Vector& b);

	// Projects a vector v1 on the line spanned by v2
	Vector Proj(const Vector& v1, const Vector& v2);

	/*Applies the Gram-Schmidt procces on the matrix A
	where columns of A are the set of vectors that needs to be made orthonormal.*/
	Matrix GramSchmidt(const Matrix& A);
}