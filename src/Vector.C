#include "Vector.H"

#include <iostream>
#include <cmath>

/* 
    Constructor
*/

Vector::Vector(double x, double y, double z) : x(x), y(y), z(z) {}

// Default Constructor
Vector::Vector() : x(0.0), y(0.0), z(0.0) {}



/* 
    Functions
*/

// Dot Product
double Vector::dotMult(const Vector &obj) const
{
    const double mult_x = x * obj.x;
    const double mult_y = y * obj.y;
    const double mult_z = z * obj.z;

    return mult_x + mult_y + mult_z;
}

// Cross Product
Vector Vector::crossMult(const Vector& obj) const
{
    const double mult_x = (y * obj.z) - (obj.y * z);
    const double mult_y = -((x * obj.z) - (obj.x * z));
    const double mult_z = (x * obj.y) - (obj.x * y);

    return Vector(mult_x, mult_y, mult_z);
}

// Modulus / Norm 
double Vector::norm() const
{
    double norm = sqrt(dotMult(*this));

    return norm;
}



/* 
    Operator Overloading
*/

// Addition +
Vector Vector::operator+(const Vector& obj) const
{
    const double z_comp = z + obj.z;
    const double x_comp = x + obj.x;
    const double y_comp = y + obj.y;

    return Vector(x_comp, y_comp, z_comp);
}

// Subtraction -
Vector Vector::operator-(const Vector& obj) const
{
    const double x_comp = x - obj.x;
    const double y_comp = y - obj.y;
    const double z_comp = z - obj.z;

    return Vector(x_comp, y_comp, z_comp);
}

// Multiplication *
Vector Vector::operator*(const double& a) const
{
    const double x_comp = x * a;
    const double y_comp = y * a;
    const double z_comp = z * a;

    return Vector(x_comp, y_comp, z_comp);
}

// Division /
Vector Vector::operator/(const double& a) const
{
    const double x_comp = x / a;
    const double y_comp = y / a;
    const double z_comp = z / a;

    return Vector(x_comp, y_comp, z_comp);
}

// Square Brackets [] 
double Vector::operator[](const int a) const
{
    if (a == 0)
    {
        return x;
    }
    else if (a == 1)
    {
        return y;
    }
    else if (a == 2)
    {
        return z;
    }
    else 
    {
        std::cerr << "Error! Incorrect indexing of Vector class object!" << std::endl;
        exit(1);
    }
}

// Not equal to !=
bool Vector::operator!=(const Vector& obj) const
{
    return x != obj.x || y != obj.y || z != obj.z;
}

// Is equal to ==
bool Vector::operator==(const Vector& obj) const
{
    return x == obj.x && y == obj.y && z == obj.z;
}