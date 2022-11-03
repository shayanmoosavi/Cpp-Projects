#include "vector.h"
#include <iostream>
#include <iomanip>

namespace LinAlg {
	// EXCEPTIONS
	VectorExceptions::VectorExceptions(const char* msg) : message(msg) {}

	const char* VectorExceptions::what() const throw() {
		return message;
	}

	// Vector Constructor
	Vector::Vector(int row) : Row(row) {
		Elem = new double[row];
	}

	// Vector Destructor (Deleted for now)
	/*Vector::~Vector() {
		delete[] Elem;
		Elem = NULL;
	}*/

	// Asking the user to enter the vector elements
	void Vector::EnterElements() {
		for (int i = 0; i < Row; i++) {
			std::cout << "Enter element " << (i + 1) << ":\n";
			std::cin >> Elem[i];
		}
	}

	// Displaying the Vector
	void Vector::Display() const {
		// Setting the precision after the decimal places
		std::cout << std::setprecision(6) << std::fixed;

		for (int i = 0; i < Row; i++) {
			std::cout << Elem[i] << std::endl;
		}
	}

	// Displaying the size of a vector
	void Vector::Size() const {
		std::cout << Row << 'x' << Col << std::endl;
	}

	// Defining the copy assignment operator
	Vector& Vector::operator=(const Vector& v) {
		this->Row = v.Row;
		for (int i = 0; i < this->Row; i++) {
			this->Elem[i] = v.Elem[i];
		}
		return *this;
	}
}

