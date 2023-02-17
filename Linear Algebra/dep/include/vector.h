#pragma once
#include <exception>

// Linear Algebra Namespace
namespace LinAlg {
    // Custom exception class for vectors
	class [[maybe_unused]] VectorExceptions : public std::exception {

        // Error Type
        //Represents what type of error caused an exception.
        const char *ErrorType;
        // Message
        // Explains the cause of the error.
        const char *Message;
	public:

        // The constructor of the VectorExceptions class
        [[maybe_unused]] explicit VectorExceptions(const char* errtyp, const char* msg);

        // Overriding the "what()" method of the base exceptions class
		[[nodiscard]] const char* what() const noexcept override;
	};

	// Vector Class
	class Vector {
	public:
		const int Col = 1; // Column Number
		int Row; // Row Number
		double* Elem; // Vector Elements


		// The constructor for the Vector class
        //Creates a Vector object with the given row as the argument.
		explicit Vector(int row);

		// The destructor for the Vector class
        //Deletes the Vector object.
		~Vector();

		// Asks the user to enter the vector elements
        [[maybe_unused]] void EnterElements() const;

		// Displays the Vector in the terminal
        [[maybe_unused]] void Display() const;

		// Displays the size of a vector
        [[maybe_unused]] void Size() const;

		// Overloading the "+" operator for vector addition
        //Adds the two corresponding vectors.
		Vector operator+(const Vector& v) const;

		// Overloading the "-" operator for vector subtraction
        //Subtracts the two corresponding vectors.
        Vector operator-(const Vector& v) const;

		// Defining the copy assignment operator
        //Copies the vector object and assigns it to the left-hand side variable.
		Vector& operator=(const Vector& v);
	};
}