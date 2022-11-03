#pragma once
#include <exception>

namespace LinAlg {
	// EXCEPTIONS
	class VectorExceptions : public std::exception {
		const char* message;
	public:

		VectorExceptions(const char* msg);

		const char* what() const throw();
	};

	// DECLARATION OF VECTOR CLASS

	// Vector Class
	class Vector {
	public:
		const int Col = 1; // Column Number
		int Row; // Row Number
		double* Elem; // Vector Elements


		// Vector Constructor
		Vector(int row);

		// Vector Destructor
		//~Vector();

		// Asking the user to enter the vector elements
		void EnterElements();

		// Displaying the Vector
		void Display() const;

		// Displaying the size of a vector
		void Size() const;

		// Overloading the "+" operator for vector addition
		Vector operator+(const Vector& v);

		// Overloading the "-" operator for vector subtraction
		Vector operator-(const Vector& v);

		// Defining the copy assignment operator
		Vector& operator=(const Vector& v);
	};
}