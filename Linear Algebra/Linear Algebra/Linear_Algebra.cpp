// Linear Algebra.cpp : Defines the entry point for the application.
//

#include "Linear Algebra.h"
#include <cmath>

// SECTION 1 : VECTOR AND MATRIX CLASSES
// ----------------------------------------------------------


// Vector Class
class Vector {
public:
	const int Col = 1; // Column Number
	int Row; // Row Number
	double* Elem; // Vector Elements


	// Vector Constructor
	Vector(int row) {
		Row = row; // Initializing the row number
		Elem = new double[row]; // Dynamically allocating the memory
	}

	// Asking the user to enter the vector elements
	void EnterElements() {
		for (int i = 0; i < Row; i++) {
			std::cout << "Enter element " << (i + 1) << ":\n";
			std::cin >> Elem[i];
		}
	}

	// Displaying the Vector
	void Display() {
		for (int i = 0; i < Row; i++) {
			std::cout << Elem[i] << std::endl;
		}
	}

	// Displaying the size of a vector
	void Size() {
		std::cout << Row << 'x' << Col << std::endl;
	}

	// Overloading the "+" operator for vector addition
	Vector operator+(Vector& v) {
		if (Row != v.Row) {
			std::logic_error vec_add_error("Undefined Operation! Vector dimensions don\'t match.");
			throw vec_add_error;
		}
		else {
			Vector UplusV(Row);
			for (int i = 0; i < Row; i++) {
				UplusV.Elem[i] = Elem[i] + v.Elem[i];
			}
			return UplusV;
		}
	}

	// Overloading the "-" operator for vector subtraction
	Vector operator-(Vector& v) {
		if (Row != v.Row) {
			std::logic_error vec_sub_error("Undefined Operation! Vector dimensions don\'t match.");
			throw vec_sub_error;
		}
		else {
			Vector UplusV(Row);
			for (int i = 0; i < Row; i++) {
				UplusV.Elem[i] = Elem[i] - v.Elem[i];
			}
			return UplusV;
		}
	}

	// Defining the copy assignment operator
	Vector& operator=(Vector& v) {
		this->Row = v.Row;
		for (int i = 0; i < this->Row; i++) {
			this->Elem[i] = v.Elem[i];
		}
		return *this;
	}
};


// Matrix Class
class Matrix {
	public:
		int Row; // Row Number
		int Col; // Column Number
		double** Elem; // Matrix Elements
				
		// Matrix Constructor
		Matrix(int row, int col) {
			
			// Initializing rows and columns
			Row = row;
			Col = col;

			// Dynamically Allocating Memory
			Elem = new double* [row];
			for (int i = 0; i < row; i++) {
				Elem[i] = new double[col];
			}
			
		}


		// Defining vector to matrix conversion
		Matrix(Vector& v) {
			
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
		void EnterElements() {
			
			// Entering Matrix Elemnts
			for (int i = 0; i < Row; i++) {
				std::cout << "Enter row " << (i + 1) << ":\n";
				for (int j = 0; j < Col; j++) {
					std::cin >> Elem[i][j];
				}
			}
		}

		// Displaying the matrix
		void Display() {
			
			// Displaying the matrix elements
			for (int i = 0; i < Row; i++) {
				for (int j = 0; j < Col; j++) {
					std::cout << Elem[i][j] << " ";
				}
				std::cout << std::endl;
			}
		}

		// Displaying the size of a matrix
		void Size() {
			std::cout << Row << 'x' << Col << std::endl;
		}

		// Overloading the "+" operator for matrix addition
		Matrix operator+(Matrix& A) {
			// Matrix dimensions should match
			if (Row != A.Row || Col != A.Col) {

				// Producing an error
				std::logic_error sum_error("Undefined Operation! Matrix dimensions don\'t match.");
				throw sum_error;
			}

			// Creating a new matrix object
			Matrix AplusB(A.Row, A.Col);

			// Doing the matrix addition
			for (int i = 0; i < A.Row; i++) {
				for (int j = 0; j < A.Col; j++) {
					AplusB.Elem[i][j] = Elem[i][j] + A.Elem[i][j];
				}
			}

			return AplusB;
		}

		// Overloading the "-" operator for matrix subtraction
		Matrix operator-(Matrix& A) {
			// Matrix dimensions should match
			if (Row != A.Row || Col != A.Col) {

				// Producing an error
				std::logic_error sub_error("Undefined Operation! Matrix dimensions don\'t match.");
				throw sub_error;
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
		double Trace() {
			
			// trace of a non-square matrix is undefined
			if (Row != Col) {
				
				// producing an error
				std::logic_error trace_error("Undefined operation! Trace is undefined for non-square matrices.");
				throw trace_error;
			}
			else {
				
				double sum = 0;
				for (int i = 0; i < Col; i++) {
					sum += Elem[i][i];
				}
				
				return sum;
			}
			
		}

		// Calsulating the determinant of a matrix
		double Det() {

			if (Row != Col) {
				std::logic_error det_error("ERROR: Determinant is undefined for non-square matrices.");
				throw det_error;
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
		Matrix Transpose() {
			
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
		Matrix SubMatrix(int row, int col) {
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
		Matrix MatAugment(Matrix& A) {
			if (Row != A.Row) {
				std::logic_error augment_error("ERROR: Rows do not match.");
				throw augment_error;
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
		void Split(int col, Matrix* Mat1, Matrix* Mat2) {
			if (col == 0) {
				throw std::logic_error("ERROR: Column index should be at least 1!");
			}
			else if (col == Col) {
				throw std::logic_error("ERROR: Can\'t split at the last column!");
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
		/*void GetCols() {

		}*/

		// Checking whether the matrix is symmetric
		bool isSymmetric() {
			
			// Symmetry is only defined for square matrices 
			if (Col != Row) {
				// producing an error
				throw std::logic_error("Symmetry is undefined for non-square matrices!");
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
		bool isSquare() {
			if (Row == Col) {
				return true;
			}
			else {
				return false;
			}
		}
};

// SECTION 2 : FUNCTIONS AND BASIC ALGORITHMS
// ----------------------------------------------------------

// Defining matrix to vector conversion
Vector ConvertToMatrix(Matrix& A) {
	if (A.Col != 1) {
		std::logic_error conversion_error("ERROR: Can\'t convert a non-column matrix to a vector.");
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

// overloading the "*" operator for matrix scaler multiplication
Matrix operator*(double c, Matrix& A) {
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

// Calculating the matrix addition
Matrix MatSum(Matrix& A, Matrix& B) {

	// Matrix dimensions should match
	if (A.Row != B.Row || A.Col != B.Col) {

		// Producing an error
		std::logic_error sum_error("Undefined Operation! Matrix dimensions don\'t match.");
		throw sum_error;
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

// Calculating the matrix multiplication
Matrix MatMult(Matrix& A, Matrix& B) {

	if (A.Col != B.Row) {

		std::logic_error mult_error("Undefined Operation! The rows of A does not match columns of B.");
		throw mult_error;
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
Matrix MatScalerMult(double c, Matrix& A) {

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

// Multiplying a row by a number
Matrix RowMult(Matrix& A, int row, double num) {
	for (int i = 0; i < A.Col; i++) {
		A.Elem[row][i] = num * A.Elem[row][i];
	}
	return A;
}

// Swapping two rows
Matrix RowSwap(Matrix& A, int row1, int row2) {
	for (int i = 0; i < A.Col; i++) {
		double temp;
		temp = A.Elem[row1][i];
		A.Elem[row1][i] = A.Elem[row2][i];
		A.Elem[row2][i] = temp;
	}
	return A;
}

// Adding row2 multiplied by a number to row1
Matrix RowAdd(Matrix& A, int row1, int row2, double num = 1) {
	for (int i = 0; i < A.Col; i++) {
		A.Elem[row1][i] += num * A.Elem[row2][i];
	}
	return A;
}

// Calculating the inverse of a matrix
Matrix Inv(Matrix& A) {
	if (A.Det() == 0) {
		// Singular Matrices are irreversible
		throw std::logic_error("ERROR: Matrix is singular!");
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
Matrix VecTranspose(Vector& v) {
	Matrix TransMat(1, v.Row);
	for (int i = 0; i < v.Row; i++) {
		TransMat.Elem[0][i] = v.Elem[i];
	}
	return TransMat;
}

// Overloading the "*" operator for vector scaler multiplication
Vector operator*(double c, Vector& v) {
	Vector cv(v.Row);
	for (int i = 0; i < v.Row; i++) {
		cv.Elem[i] = c * v.Elem[i];
	}
	return cv;
}

// Adding two vectors
Vector Vecadd(Vector& u, Vector& v) {
	if (u.Row != v.Row) {
		std::logic_error vec_add_error("Undefined Operation! Vector dimensions don\'t match.");
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

// Multiplying a vector by a scaler
Vector VecScalerMult(double c, Vector& v) {
	Vector cv(v.Row);
	for (int i = 0; i < v.Row; i++) {
		cv.Elem[i] = c * v.Elem[i];
	}
	return cv;
}

// Calculating the cross product
Vector Cross(Vector& u, Vector& v) {
	if (u.Row != v.Row) {
		throw std::logic_error("Undefined Operation! Vector dimensions don\'t match.");
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
double Dot(Vector& u, Vector& v) {
	Matrix V = v;
	return MatMult(VecTranspose(u), V).Elem[0][0];
}

/*Solves the linear system Ax=b where A is a matrix, 
x is an unknown vector and b is a vector*/
Matrix LinearSystemSolve(Matrix& A, Vector& b) {
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
Vector Proj(Vector& v1, Vector& v2) {
	return VecScalerMult(Dot(v1, v2) / pow(Dot(v2, v2), 0.5), v2);
}

/*Applies the Gram-Schmidt procces on the matrix A
where matrix A is the set of vectors that needs to be made orthonormal.*/ 
Matrix GramSchmidt(Matrix& A) {
	
	// Initializing the vectors and a matrix to store our answer in
	Vector* v1 = new Vector(A.Row);
	Vector* v2 = new Vector(A.Row);
	Vector* vPerp = new Vector(A.Row);
	Matrix B(A.Row, A.Col);
	
	// GRAM-SCHMIDT PROCESS ALGORITHM
	// Setting the vector elements from the input matrix 
	for (int i = 0; i < A.Col; i++) {
		for (int j = 0; j < A.Row; j++) {
			v1->Elem[j] = A.Elem[j][i];
		}
		// Debugging code:
		v1->Display();
		if (i == 0) {
			*v1 = (1 / pow(Dot(*v1, *v1), 0.5)) * (*v1); // Normalizing the vector
			*v2 = *v1;
			for (int j = 0; j < A.Row; j++) {
				B.Elem[j][i] = v1->Elem[j];
			}
			std::cout << "vector " << i + 1 << std::endl;
			v2->Display();
		}
		else {
			// Making an orthogonal vector to the previous vector
			*vPerp = *v1;
			for (int j = i; j > 0; j--) {
				for (int k = 0; k < A.Row; k++) {
					v2->Elem[k] = B.Elem[k][j];
				}
				*vPerp = *v1 - Proj(*v1, *v2);
				*v2 = *vPerp;
				*v2 = (1 / pow(Dot(*v2, *v2), 0.5)) * (*v2);
				for (int k = 0; k < A.Row; k++) {
					B.Elem[k][j] = v2->Elem[k];
				}
				std::cout << "vector " << i + 1 << std::endl;
				v2->Display();
			}
			/*if (i == 1) {
				*vPerp = *v1 - Proj(*v1, *v2);
				*v2 = *vPerp;
				*v2 = (1 / pow(Dot(*v2, *v2), 0.5)) * (*v2);
				for (int j = 0; j < A.Row; j++) {
					B.Elem[j][i] = v2->Elem[j];
					std::cout << "vector " << i + 1 << std::endl;
					v2->Display();
				}
			}
			else {
				
			}*/
			
		}
	}
	
	return B;
}

int main()
{
	//Matrix A(2, 3); // Creating a 2x2 matrix
	//A.EnterElements(); // Entering matrix elemnents
	//std::cout << "Before transpose:\n";
	//std::cout << "Size: ";
	//A.Size();
	//A.Display();
	//std::cout << "After transpose:\n";
	//std::cout << "Size: ";
	//Matrix B = A.Transpose();
	//B.Size();
	//B.Display();
	//
	//// Deallocating the memory
	//delete A.Elem;
	//A.Elem = NULL;
	//delete B.Elem;
	//B.Elem = NULL;

	/*Matrix A(3, 3);
	int col = 2;
	Matrix* A1 = new Matrix(3, col);
	Matrix* A2 = new Matrix(3, A.Col - col);
	A.EnterElements();
	A.Display();
	A.Split(col, A1, A2);
	A1->Display();
	A2->Display();

	delete A1;
	delete A2;
	A1 = NULL;
	A2 = NULL;*/

	/*Matrix A(3, 3);
	A.EnterElements();
	A.Display();
	Inv(A).Display();*/

	/*Vector v(3);
	v.EnterElements();
	std::cout << typeid(v).name();
	Matrix A(3, 1);
	A = v;
	A.Display();
	std::cout << typeid(A).name();*/
	
	/*Vector v(3);
	v.EnterElements();
	Matrix A = VecTranspose(v);
	A.Display();*/

	/*Vector b(3);
	b.EnterElements();
	b.Display();
	Matrix A(3, 3);
	A.EnterElements();
	A.Display();
	Matrix x = LinearSystemSolve(A, b);
	std::cout << "The solution to the system is: \n";
	x.Display();*/

	/*Vector v1(3), v2(3);
	v1.EnterElements();
	v2.EnterElements();
	Proj(v1, v2).Display();*/

	/*Matrix A(2, 2), B(2, 2);
	A.EnterElements();
	B.EnterElements();
	Matrix C = A + B;
	C.Display();*/

	Matrix A(3, 3);
	A.EnterElements();
	GramSchmidt(A).Display();


	return 0;
}
