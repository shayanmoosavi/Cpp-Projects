#include "linear_algebra.h"
#include <cmath>


// Linear Algebra namespace
namespace LinAlg {
    [[maybe_unused]] Vector ConvertToMatrix(const Matrix& A) {
		if (A.Col != 1) {
			throw VectorExceptions("[ConversionException]",
                                   "ERROR: Cannot convert a non-column matrix to a vector.");
		}
		else {
			static Vector v(A.Row);
			for (int i = 0; i < v.Row; i++) {
				v.Elem[i] = A.Elem[i][0];
			}
			return v;
		}
	}

	Matrix Identity(int d) {

        static Matrix I(d, d);

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

    [[maybe_unused]] Matrix Diagonal(const int d, const double var[]) {

		static Matrix DiagMat(d, d);

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

    [[maybe_unused]] Matrix MatSum(const Matrix& A, const Matrix& B) {

		// Matrix dimensions should match
		if (A.Row != B.Row || A.Col != B.Col) {

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
				AplusB.Elem[i][j] = A.Elem[i][j] + B.Elem[i][j];
			}
		}

		return AplusB;
	}

	Matrix MatMult(const Matrix& A, const Matrix& B) {

		if (A.Col != B.Row) {
			throw MatrixExceptions("[NonMatchingDimensionsException]",
                                   "Undefined Operation! The rows of A does not match columns of B.");
		}
		else {
            /* Creating a new matrix object. Declared static in order to prevent the destructor
         from being called inappropriately.*/
			static Matrix AB(A.Row, B.Col);

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

	Matrix MatScalerMult(double c, const Matrix& A) {

        /* Creating a new matrix object. Declared static in order to prevent the destructor
         from being called inappropriately.*/
		static Matrix cA(A.Row, A.Col);

		// Doing the scaler multiplication
		for (int i = 0; i < A.Row; i++) {
			for (int j = 0; j < A.Col; j++) {
				cA.Elem[i][j] = c * A.Elem[i][j];
			}
		}

		return cA;
	}

	Matrix operator*(double c, const Matrix& A) {

		return MatScalerMult(c, A);
	}

    [[maybe_unused]] Matrix RowMult(const Matrix& A, int row, double num) {
		for (int i = 0; i < A.Col; i++) {
			A.Elem[row][i] = num * A.Elem[row][i];
		}

		return A;
	}

    [[maybe_unused]] Matrix RowSwap(const Matrix& A, int row1, int row2) {
		for (int i = 0; i < A.Col; i++) {
			double temp;
			temp = A.Elem[row1][i];
			A.Elem[row1][i] = A.Elem[row2][i];
			A.Elem[row2][i] = temp;
		}

		return A;
	}

    [[maybe_unused]] Matrix RowAdd(const Matrix& A, int row1, int row2, double num = 1) {
		for (int i = 0; i < A.Col; i++) {
			A.Elem[row1][i] += num * A.Elem[row2][i];
		}

		return A;
	}

    [[maybe_unused]] Matrix Inv(Matrix& A) {
		if (A.Det() == 0) {
			// Singular Matrices are irreversible
			throw MatrixExceptions("[SingularMatrixException]",
                                   "ERROR: Matrix is singular!");

		}
		else {

			static Matrix B = A.MatAugment(Identity(A.Row));
			auto* B1 = new Matrix(A.Row, A.Col);
			auto* B2 = new Matrix(A.Row, A.Col);

            // GAUSS JORDAN ELIMINATION ALGORITHM
            for (int i = 0; i < A.Row; i++) {

                // Swapping the row having diagonal element of zero with another row
                if (B.Elem[i][i] == 0) {

                    if (i != A.Col - 1) {

                        for (int j = 0; j < B.Col; j++) {

                            double temp;
                            temp = B.Elem[i][j];
                            B.Elem[i][j] = B.Elem[i+1][j];
                            B.Elem[i+1][j] = temp;
                        }
                    }
                    else {

                        // Last row will exceed array bounds, so we swap with the previous row
                        for (int j = 0; j < B.Col; j++) {

                            double temp;
                            temp = B.Elem[i][j];
                            B.Elem[i][j] = B.Elem[i-1][i];
                            B.Elem[i-1][j] = temp;
                        }
                    }
                }
                else {

                    if (B.Elem[i][i] != 1) { // Skipping the row that already has a one in its pivot

                        double temp = B.Elem[i][i];

                        for (int j = 0; j < B.Col; j++) {
                            B.Elem[i][j] = 1 / temp * B.Elem[i][j];
                        }
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
                            double temp = -B.Elem[j][i];
                            for (int k = 0; k < B.Col; k++) {
                                B.Elem[j][k] += temp * B.Elem[i][k];
                            }
                        }
                    }
                }
            }

            // Splitting the augmented matrix and recovering the inverse matrix
            B.Split(A.Col, B1, B2);

            // Matrix on the left is not needed anymore so the memory can be freed
            delete B1;

            return *B2;
        }
	}

	Matrix VecTranspose(const Vector& v) {

		static Matrix TransMat(1, v.Row);

		for (int i = 0; i < v.Row; i++) {
			TransMat.Elem[0][i] = v.Elem[i];
		}

		return TransMat;
	}

	Vector Vecadd(const Vector& u, const Vector& v) {

		if (u.Row != v.Row) {
			throw VectorExceptions("[UndefinedOperationException]",
                                   "Undefined Operation! Vector dimensions do not match.");
		}
		else {

			static Vector UplusV(u.Row);

			for (int i = 0; i < u.Row; i++) {
				UplusV.Elem[i] = u.Elem[i] + v.Elem[i];
			}

			return UplusV;
		}
	}

	// Subtracting two vectors
	Vector Vecsub(const Vector& u, const Vector& v) {

		if (u.Row != v.Row) {
            throw VectorExceptions("[UndefinedOperationException]",
                                   "Undefined Operation! Vector dimensions do not match.");
        }
		else {

			static Vector UminusV(u.Row);

			for (int i = 0; i < u.Row; i++) {
				UminusV.Elem[i] = u.Elem[i] - v.Elem[i];
			}

			return UminusV;
		}
	}

	Vector Vector::operator+(const Vector& v) const {

		return Vecadd(*this, v);
	}

	Vector Vector::operator-(const Vector& v) const {

		return Vecsub(*this, v);
	}

	Vector VecScalerMult(double c, const Vector& v) {

        static Vector cv(v.Row);

        for (int i = 0; i < v.Row; i++) {
			cv.Elem[i] = c * v.Elem[i];
		}

		return cv;
	}

	Vector operator*(double c, const Vector& v) {
		return VecScalerMult(c, v);
	}

    [[maybe_unused]] Vector Cross(const Vector& u, const Vector& v) {

		if (u.Row != v.Row) {
            throw VectorExceptions("[UndefinedOperationException]",
                                   "Undefined Operation! Vector dimensions do not match.");
        }
		else if (u.Row != 3) {
            throw VectorExceptions("[UndefinedOperationException]",
                                   "Undefined Operation! Cross product is undefined for non-3D vectors.");
        }
		else {

			static Vector UcrossV(3);

			UcrossV.Elem[0] = u.Elem[1] * v.Elem[2] - u.Elem[2] * v.Elem[1];
			UcrossV.Elem[1] = u.Elem[2] * v.Elem[0] - u.Elem[0] * v.Elem[2];
			UcrossV.Elem[2] = u.Elem[0] * v.Elem[1] - u.Elem[1] * v.Elem[0];

			return UcrossV;
		}
	}

	double Dot(const Vector& u, const Vector& v) {
		auto V = Matrix(v);
		return MatMult(VecTranspose(u), V).Elem[0][0];
	}

	/*Solves the linear system Ax=b where A is a matrix,
	x is an unknown vector and b is a vector*/
    [[maybe_unused]] Matrix LinearSystemSolve(Matrix& A, const Vector& b) {
		if (A.Det() == 0) {
			// Singular Matrices are irreversible
            throw MatrixExceptions("[SingularMatrixException]",
                                   "ERROR: Matrix is singular!");
        }
		else {

			static Matrix bMat = (const Matrix) b;
			static Matrix B = A.MatAugment(bMat);
			auto* B1 = new Matrix(B.Row, B.Row);
			auto* B2 = new Matrix(b.Row, b.Col);

            // GAUSS JORDAN ELIMINATION ALGORITHM
            for (int i = 0; i < A.Row; i++) {

                // Swapping the row having diagonal element of zero with another row
                if (B.Elem[i][i] == 0) {

                    if (i != A.Col - 1) {

                        for (int j = 0; j < B.Col; j++) {

                            double temp;
                            temp = B.Elem[i][j];
                            B.Elem[i][j] = B.Elem[i+1][j];
                            B.Elem[i+1][j] = temp;
                        }
                    }
                    else {

                        // Last row will exceed array bounds, so we swap with the previous row
                        for (int j = 0; j < B.Col; j++) {

                            double temp;
                            temp = B.Elem[i][j];
                            B.Elem[i][j] = B.Elem[i-1][i];
                            B.Elem[i-1][j] = temp;
                        }
                    }
                }
                else {
                    if (B.Elem[i][i] != 1) { // Skipping the row that already has a one in its pivot

                        double temp = B.Elem[i][i];
                        for (int j = 0; j < B.Col; j++) {
                            B.Elem[i][j] = 1 / temp * B.Elem[i][j];
                        }
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

                            double temp = -B.Elem[j][i];
                            for (int k = 0; k < B.Col; k++) {
                                B.Elem[j][k] += temp * B.Elem[i][k];
                            }
                        }
                    }
                }
            }

            // Splitting the augmented matrix and recovering the inverse matrix
            B.Split(A.Col, B1, B2);

            // Matrix on the left is not needed anymore so the memory can be freed
            delete B1;

            return *B2;
        }
	}

	Vector Proj(const Vector& v1, const Vector& v2) {
		return (Dot(v1, v2) / pow(Dot(v2, v2), 0.5)) * v2;
	}

    [[maybe_unused]] Matrix GramSchmidt(const Matrix& A) {

		/* Initializing the vectors and a matrix to store our answer in.
		We use pointers in order to manipulate the vector and matrix elements directly*/
		auto* v1 = new Vector(A.Row); // New Vector Pointer
		auto* v2 = new Vector(A.Row); // Old Vector Pointer
		//Vector* vPerp = new Vector(A.Row); // Orthogonal Vector Pointer
		auto* B = new Matrix(A.Row, A.Col); // Matrix pointer to store the answer

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
				}

				*v2 = (1 / pow(Dot(*v2, *v2), 0.5)) * (*v2);

				for (int k = 0; k < A.Row; k++) {
					B->Elem[k][i] = v2->Elem[k];
				}
			}
		}
		// Freeing the memory
		delete v1;
		delete v2;

        return *B;
	}
}
