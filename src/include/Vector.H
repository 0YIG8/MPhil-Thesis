#ifndef VECTOR_H
#define VECTOR_H

class Vector
{
    private:
        double x, y, z;

    public:
        Vector(double x, double y, double z); // Constructor
        Vector(); // Default constructor

        double dotMult(const Vector& obj) const; 
        Vector crossMult(const Vector& obj) const;
        double norm() const;

        Vector operator+(const Vector& obj) const;
        Vector operator-(const Vector& obj) const;

        Vector operator*(const double& a) const; // Mult vector by scalar
        Vector operator/(const double& a) const; // Div vector by scalar
        double operator[](const int a) const; // v1[0] returns x component, v1[1] returns y, v1[2] returns z
        bool operator!=(const Vector& obj) const;
        bool operator==(const Vector& obj) const;
};

#endif