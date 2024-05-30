#include "Point.H"
#include <fstream>
#include <sstream>
#include <iostream>


int Point::nPoints = 0; // Initialize no of pts to 0

/* 
    Constructor
*/

Point::Point(double x, double y, double z)
{
    this -> x = x;
    this -> y = y;
    this -> z = z;

    index = nPoints;
    nPoints++;
}



/* 
    Functions
*/

// Convert Point to Vector type
Vector Point::pointToVector() const
{
    Vector point(x, y, z);

    return point;
}



/* 
    Get Functions
*/

// Return point from index
Point Point::getPointFromIndex(int idx, const std::vector<Point> listOfAllPoints)
{
    Point p = listOfAllPoints[idx];
    return p;
}

int Point::getPointIndex() const
{
    return index;
}



/*
    OPERATOR OVERLOADING
*/

// Not equal to !=
bool Point::operator!=(const Point& obj) const
{
    return x != obj.x || y != obj.y || z != obj.z;
}

// Is equal to ==
bool Point::operator==(const Point& obj) const
{
    return x == obj.x && y == obj.y && z == obj.z;
}
