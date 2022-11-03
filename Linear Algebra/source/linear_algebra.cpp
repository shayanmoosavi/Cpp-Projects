#include "linear_algebra.h"
#include <iostream>
#include <cmath>

namespace LinAlg {
	Vector ConvertToMatrix(const Matrix& A) {
		if (A.Col != 1) {
			std::logic_error conversion_error("ERROR: Cannot convert a non-column matrix to a vector.");
			throw conversion_error;
		}
		else {
			Vector v(A.Row);
			for (int i = 0; i < v.Row; i++) {
				v.Elem[i] = A.Elem[i][0];
			}
			return v;
		}
	}

	// Creating an identity matrix with size dxd
	Matrix Identity(int d) {
		Matrix I(d, d);
		for (int i = 0; i < d; i++) {
			for (int j = 0; j < d; j++) {
				if (i == j) {
					I.Elem[i][j] = 1;
				}
				else {
					I.Elem[i][j] = 0;
				}
			}
		}

		return I;
	}

	// Creating a diagonal matrix with diagonal elements
	Matrix Diagonal(const int d, double var[]) {
		Matrix DiagMat(d, d);
		var[d];
		for (int i = 0; i < d; i++) {
			for (int j = 0; j < d; j++) {
				if (i == j) {
					DiagMat.Elem[i][j] = var[i];
				}
				else {
					DiagMat.Elem[i][j] = 0;
				}
			}
		}
		return DiagMat;
	}

	// Calculating the matrix addition
	Matrix MatSum(const Matrix& A, const Matrix& B) {

		// Matrix dimensions should match
		if (A.Row != B.Row || A.Col != B.Col) {

			// Producing an error
			MatrixExceptions sum_err("Undefined Operation! Matrix dimensions do not match.");
			throw sum_err;
		}

		// Creating a new matrix object
		Matrix AplusB(A.Row, A.Col);

		// Doing the matrix addition
		for (int i = 0; i < A.Row; i++) {
			for (int j = 0; j < A.Col; j++) {
				AplusB.Elem[i][j] = A.Elem[i][j] + B.Elem[i][j];
			}
		}

		return AplusB;
	}

	// Overloading the "+" operator for matrix addition
	Matrix Matrix::operator+(const Matrix& A) {

		return MatSum(*this, A);
	}


	// Calculating the matrix multiplication
	Matrix MatMult(const Matrix& A, const Matrix& B) {

		if (A.Col != B.Row) {

			MatrixExceptions mult_err("Undefined Operation! The rows of A does not match columns of B.");
			throw mult_err;
		}
		else {

			// Creating a new matrix object
			Matrix AB(A.Row, B.Col);

			// Initializing the matrix elements to zero
			for (int i = 0; i < AB.Row; i++) {
				for (int j = 0; j < AB.Col; j++) {
					AB.Elem[i][j] = 0;
				}
			}

			// Doing the matrix multiplication
			for (int i = 0; i < AB.Row; i++) {
				for (int j = 0; j < AB.Col; j++) {
					for (int k = 0; k < A.Col; k++) {
						AB.Elem[i][j] += A.Elem[i][k] * B.Elem[k][j];
					}
				}
			}

			return AB;
		}
	}

	// Calculating the scaler c multiplied by natrix A
	Matrix MatScalerMult(double c, const Matrix& A) {

		// Creating a new matrix object
		Matrix cA(A.Row, A.Col);

		// Doing the scaler multiplication
		for (int i = 0; i < A.Row; i++) {
			for (int j = 0; j < A.Col; j++) {
				cA.Elem[i][j] = c * A.Elem[i][j];
			}
		}

		return cA;
	}

	// overloading the "*" operator for matrix scaler multiplication
	Matrix operator*(double c, const Matrix& A) {

		return MatScalerMult(c, A);
	}


	// Multiplying a row by a number
	Matrix RowMult(const Matrix& A, int row, double num) {
		for (int i = 0; i < A.Col; i++) {
			A.Elem[row][i] = num * A.Elem[row][i];
		}
		return A;
	}

	// Swapping two rows
	Matrix RowSwap(const Matrix& A, int row1, int row2) {
		for (int i = 0; i < A.Col; i++) {
			double temp;
			temp = A.Elem[row1][i];
			A.Elem[row1][i] = A.Elem[row2][i];
			A.Elem[row2][i] = temp;
		}
		return A;
	}

	// Adding row2 multiplied by a number to row1
	Matrix RowAdd(const Matrix& A, int row1, int row2, double num = 1) {
		for (int i = 0; i < A.Col; i++) {
			A.Elem[row1][i] += num * A.Elem[row2][i];
		}
		return A;
	}

	// Calculating the inverse of a matrix
	Matrix Inv(Matrix& A) {
		if (A.Det() == 0) {
			// Singular Matrices are irreversible
			MatrixExceptions singular_err("ERROR: Matrix is singular!");
			throw singular_err;
		}
		else {
			Matrix B = A.MatAugment(Identity(A.Row));
			Matrix* B1 = new Matrix(A.Row, A.Col);
			Matrix* B2 = new Matrix(A.Row, A.Col);

			// GAUSS JORDAN ELIMINATION ALGORITHM
			for (int i = 0; i < A.Row; i++) {
				// Swapping the row having diagonal element of zero with another row
				if (B.Elem[i][i] == 0) {
					if (i != A.Col - 1) {
						RowSwap(B, i, i + 1);
					}
					else {
						// Last row will exceed array bounds so we swap with the previous row
						RowSwap(B, i, i - 1);
					}
				}
				else {
					if (B.Elem[i][i] != 1) {
						// Skipping the row that already has a one in its pivot
						RowMult(B, i, 1 / B.Elem[i][i]);
					}
					for (int j = 0; j < B.Row; j++) {
						if (j == i) {
							continue;
						}
						else if (B.Elem[j][i] == 0) {
							// Skipping the row that already has zero under its pivot
							continue;
						}
						else {
							RowAdd(B, j, i, -B.Elem[j][i]);
						}
					}
				}
			}
			// Splitting the augmented matrix and recovering the inverse matrix
			B.Split(A.Col, B1, B2);
			// Matrix on the left is not needed anymore so the memory can be freed
			delete B1;
			B1 = NULL;
			return *B2;
		}
	}

	// Transposing a vector
	Matrix VecTranspose(const Vector& v) {
		Matrix TransMat(1, v.Row);
		for (int i = 0; i < v.Row; i++) {
			TransMat.Elem[0][i] = v.Elem[i];
		}
		return TransMat;
	}

	// Adding two vectors
	Vector Vecadd(const Vector& u, const Vector& v) {
		if (u.Row != v.Row) {
			std::logic_error vec_add_error("Undefined Operation! Vector dimensions do not match.");
			throw vec_add_error;
		}
		else {
			Vector UplusV(u.Row);
			for (int i = 0; i < u.Row; i++) {
				UplusV.Elem[i] = u.Elem[i] + v.Elem[i];
			}
			return UplusV;
		}
	}

	// Subtracing two vectors
	Vector Vecsub(const Vector& u, const Vector& v) {
		if (u.Row != v.Row) {
			std::logic_error vec_sub_error("Undefined Operation! Vector dimensions do not match.");
			throw vec_sub_error;
		}
		else {
			Vector UminusV(u.Row);
			for (int i = 0; i < u.Row; i++) {
				UminusV.Elem[i] = u.Elem[i] - v.Elem[i];
			}
			return UminusV;
		}
	}

	// Overloading the "+" operator for vector addition
	Vector Vector::operator+(const Vector& v) {

		return Vecadd(*this, v);
	}

	// Overloading the "-" operator for vector subtraction
	Vector Vector::operator-(const Vector& v) {

		return Vecsub(*this, v);
	}

	// Multiplying a vector by a scaler
	Vector VecScalerMult(double c, const Vector& v) {
		Vector cv(v.Row);
		for (int i = 0; i < v.Row; i++) {
			cv.Elem[i] = c * v.Elem[i];
		}
		return cv;
	}

	// Overloading the "*" operator for vector scaler multiplication
	Vector operator*(double c, const Vector& v) {
		return VecScalerMult(c, v);
	}

	// Calculating the cross product
	Vector Cross(const Vector& u, const Vector& v) {
		if (u.Row != v.Row) {
			throw std::logic_error("Undefined Operation! Vector dimensions do not match.");
		}
		else if (u.Row != 3) {
			throw std::logic_error("ERROR: Cross product is undefined for non-3D vectors.");
		}
		else {

			Vector UcrossV(3);
			UcrossV.Elem[0] = u.Elem[1] * v.Elem[2] - u.Elem[2] * v.Elem[1];
			UcrossV.Elem[1] = u.Elem[2] * v.Elem[0] - u.Elem[0] * v.Elem[2];
			UcrossV.Elem[2] = u.Elem[0] * v.Elem[1] - u.Elem[1] * v.Elem[0];

			return UcrossV;
		}
	}

	// Calculating the dot product
	double Dot(const Vector& u, const Vector& v) {
		Matrix V = v;
		return MatMult(VecTranspose(u), V).Elem[0][0];
	}

	/*Solves the linear system Ax=b where A is a matrix,
	x is an unknown vector and b is a vector*/
	Matrix LinearSystemSolve(Matrix& A, const Vector& b) {
		if (A.Det() == 0) {
			// Singular Matrices are irreversible
			throw std::logic_error("ERROR: Matrix is singular!");
		}
		else {
			Matrix b = b;
			Matrix B = A.MatAugment(b);
			Matrix* B1 = new Matrix(A.Row, A.Col);
			Matrix* B2 = new Matrix(b.Row, b.Col);

			// GAUSS JORDAN ELIMINATION ALGORITHM
			for (int i = 0; i < A.Row; i++) {
				// Swapping the row having diagonal element of zero with another row
				if (B.Elem[i][i] == 0) {
					if (i != A.Col - 1) {
						RowSwap(B, i, i + 1);
					}
					else {
						// Last row will exceed array bounds so we swap with the previous row
						RowSwap(B, i, i - 1);
					}
				}
				else {
					if (B.Elem[i][i] != 1) {
						// Skipping the row that already has a one in its pivot
						RowMult(B, i, 1 / B.Elem[i][i]);
					}
					for (int j = 0; j < B.Row; j++) {
						if (j == i) {
							continue;
						}
						else if (B.Elem[j][i] == 0) {
							// Skipping the row that already has zero under its pivot
							continue;
						}
						else {
							RowAdd(B, j, i, -B.Elem[j][i]);
						}
					}
				}
			}
			// Splitting the augmented matrix and recovering the inverse matrix
			B.Split(A.Col, B1, B2);
			// Matrix on the left is not needed anymore so the memory can be freed
			delete B1;
			B1 = NULL;

			return *B2;
		}
	}

	// Projects a vector v1 on the line spanned by v2
	Vector Proj(const Vector& v1, const Vector& v2) {
		return (Dot(v1, v2) / pow(Dot(v2, v2), 0.5)) * v2;
	}

	/*Applies the Gram-Schmidt procces on the matrix A
	where columns of A are the set of vectors that needs to be made orthonormal.*/
	Matrix GramSchmidt(const Matrix& A) {

		/* Initializing the vectors and a matrix to store our answer in.
		We use pointers in order to manipulate the vector and matrix elements directly*/
		Vector* v1 = new Vector(A.Row); // New Vector Pointer
		Vector* v2 = new Vector(A.Row); // Old Vector Pointer
		//Vector* vPerp = new Vector(A.Row); // Orthogonal Vector Pointer
		Matrix* B = new Matrix(A.Row, A.Col); // Matrix pointer to store the answer

		// GRAM-SCHMIDT PROCESS ALGORITHM
		// Setting the vector elements from the input matrix 
		for (int i = 0; i < A.Col; i++) {
			for (int j = 0; j < A.Row; j++) {
				v1->Elem[j] = A.Elem[j][i];
			}
			if (i == 0) {
				*v1 = (1 / pow(Dot(*v1, *v1), 0.5)) * (*v1); // Normalizing the vector
				*v2 = *v1;
				for (int j = 0; j < A.Row; j++) {
					B->Elem[j][i] = v1->Elem[j];
				}
			}
			else {
				// Making an orthogonal vector to the previous vector
				int j = i;
				while (j != 0) {
					for (int k = 0; k < A.Row; k++) {
						v2->Elem[k] = B->Elem[k][j - 1];
					}
					*v1 = *v1 - Proj(*v1, *v2);
					*v2 = *v1;
					j--;
				};
				*v2 = (1 / pow(Dot(*v2, *v2), 0.5)) * (*v2);
				for (int k = 0; k < A.Row; k++) {
					B->Elem[k][i] = v2->Elem[k];
				}
			}
		}
		// Freeing the memory
		delete v1;
		delete v2;
		v1 = NULL;
		v2 = NULL;

		return *B;
	}
}
