#include "vector.h"
#include <iostream>
#include <iomanip>
#include <cstring>


// Linear Algebra Namespace
namespace LinAlg {

    [[maybe_unused]] VectorExceptions::VectorExceptions(
            // Member Initializer List
            const char* errtyp,const char* msg) : ErrorType(errtyp), Message(msg) {}

	const char* VectorExceptions::what() const noexcept {
        // Concatenating the strings
        static char buffer[sizeof(this->ErrorType) + sizeof(" - ") + sizeof(this->Message)];
        strcpy(buffer, ErrorType);
        strcat(buffer, " - ");
        strcat(buffer, Message);

        return buffer;
	}

	Vector::Vector(int row) : Row(row) {
		// Dynamically allocating the memory
        Elem = new double[row];
	}

	Vector::~Vector() {
		delete[] Elem;
	}

    [[maybe_unused]] void Vector::EnterElements() const {
		for (int i = 0; i < Row; i++) {
			std::cout << "Enter element " << (i + 1) << ":\n";
			std::cin >> Elem[i];
		}
	}

    [[maybe_unused]] void Vector::Display() const {
		// Setting the precision after the decimal places
		std::cout << std::setprecision(6) << std::fixed;

		for (int i = 0; i < Row; i++) {
			std::cout << Elem[i] << std::endl;
		}
	}

    [[maybe_unused]] void Vector::Size() const {
		std::cout << Row << 'x' << Col << std::endl;
	}

	Vector& Vector::operator=(const Vector& v) {

        if(this != &v){ // Checking for self assignment

            this->Row = v.Row;

            for (int i = 0; i < this->Row; i++) {
                // Assigning the values to the left-hand side variables
                this->Elem[i] = v.Elem[i];
            }
        }
		return *this;
	}
}

