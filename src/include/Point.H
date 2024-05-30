#ifndef POINT_H
#define POINT_H

#include <vector>
#include <string>
#include "Vector.H"

class Point
{
    static int nPoints;  // Number of Points

    private:
        double x, y, z;
        int index;

    public:
        Point(double, double, double); // Constructor

        Vector pointToVector() const;
        
        int getPointIndex() const;
        static Point getPointFromIndex(int idx, const std::vector<Point> listOfAllPoints);

        bool operator!=(const Point& obj) const; // Overloaded !=
        bool operator==(const Point& obj) const; // Overloaded ==

};

#endif