#include "Face.H"
#include <iostream>


int Face::nFaces = 0; // Initialize no of faces to 0

/*
    Constructor
*/

Face::Face(std::vector<Point> lp)
{    
    listOfPoints = lp;

    faceArea = 0.0;
    index = nFaces;
    ownerIndex = -1;
    neighbourIndex = -1;
    delta = 0.0;
    interpolationFactor = 0.0;
    faceFlux = 0.0;

    calcFaceCentre();
    calcFaceNormal();
    calcFaceArea();

    isBoundaryFace = false;
    boundaryName = "none";

    isMovingFace = false;
    faceZone = "default";

    nFaces++;
}



/*
    Calculations
*/

// Calculate Face Centre
void Face::calcFaceCentre()
{
    Vector centroid; // Vector to store centroid

    // Calculate centroid (naive)
    for (size_t i = 0; i < listOfPoints.size(); i++)
    {
        Point temp = listOfPoints[i];
        Vector v = temp.pointToVector();
        centroid = centroid + v;
    }

    centroid = centroid / listOfPoints.size(); // Average of sum

    if (listOfPoints.size() == 3) // If the face contains only 3 points
    {
        faceCentre = centroid;
    }
    else 
    {
        // Calc area of individual triangles
        for (size_t i = 0; i < listOfPoints.size(); i++)
        {
            // Get 2 consecutive points
            Point p1 = listOfPoints[i];
            Point p2 = listOfPoints[(i + 1) % listOfPoints.size()];

            // Convert to Vector
            Vector v1 = p1.pointToVector();
            Vector v2 = p2.pointToVector();

            // Calc centroid of small triangle
            Vector triCentroid = (v1 + v2 + centroid) / 3;

            // face_centre = face_centre + (tri_centroid * tri_area);
            faceCentre = faceCentre + (triCentroid / listOfPoints.size());
        }
    }
}


// Calculate Face Area
void Face::calcFaceArea()
{
    // Calc area of individual triangles
    for (size_t i = 0; i < listOfPoints.size(); i++)
    {
        // Get 2 consecutive points
        Point p1 = listOfPoints[i];
        Point p2 = listOfPoints[(i + 1) % listOfPoints.size()];

        // Convert to Vector
        Vector v1 = p1.pointToVector();
        Vector v2 = p2.pointToVector();

        // Get vectors to centroid and calc cross product to get area of //gm
        Vector cross = (v1 - faceCentre).crossMult(v2 - faceCentre);

        // Add ((norm of //gm area) / 2) to total face area
        faceArea = faceArea + (cross.norm() * 0.5);
    }
}


// Calculate face normal
void Face::calcFaceNormal()
{
    // Calc area of individual triangles
    for (size_t i = 0; i < listOfPoints.size(); i++)
    {
        // Get 2 consecutive points
        Point p1 = listOfPoints[i];
        Point p2 = listOfPoints[(i + 1) % listOfPoints.size()];

        // Convert to Vector
        Vector v1 = p1.pointToVector();
        Vector v2 = p2.pointToVector();

        // Get vectors to centroid and cross mult to get //gm area, then divide by 2
        Vector cross = ((v1 - faceCentre).crossMult(v2 - faceCentre)) * 0.5;
        
        faceNormal = faceNormal + cross;
    }
    
    faceNormal = faceNormal / faceNormal.norm();
}



/* 
    Set Functions
*/

// Set face owner
void Face::setFaceOwner(int index)
{
    ownerIndex = index;
}

// Set face neighbour
void Face::setFaceNeighbour(int index)
{
    neighbourIndex = index;
}

// Set a face as a boundary face
void Face::setAsBoundaryFace()
{
    isBoundaryFace = true; 
}

// Set name of boundary face
void Face::setBoundaryName(std::string boundaryName)
{
    this->boundaryName = boundaryName;
}

// Set face delta
void Face::setFaceDelta(double deltaValue)
{
    delta = deltaValue;
}

// Set face interpolation factor
void Face::setFaceInterpolationFactor(double intFactValue)
{
    interpolationFactor = intFactValue;
}

// Set face flux
void Face::setFaceFlux(double flux)
{
    faceFlux = flux;
}

// Set face as a moving face
void Face::setAsMovingFace()
{
    isMovingFace = true; 
}

// Set face zone name
void Face::setFaceZone(std::string zoneName)
{
    faceZone = zoneName;
}



/*
    Get Functions
*/

// Get Face Centre
Vector Face::getFaceCentre() const
{
    return faceCentre;
}

// Get Face Area
double Face::getFaceArea() const
{
    return faceArea;
}

// Get Face Normal
Vector Face::getFaceNormal() const
{
    return faceNormal;
}

// Get List of Points
std::vector<Point> Face::getListOfPoints() const
{
    return listOfPoints;
}

// Get owner 
int Face::getFaceOwner() const
{
    return ownerIndex;
}

// Get neighbour 
int Face::getFaceNeighbour() const
{
    return neighbourIndex;
}

// Get face index
int Face::getFaceIndex() const
{
    return index;
}

// Return true if a face is a boundary face
bool Face::getIfBoundaryFace() const
{
    return isBoundaryFace;
}

// Get face boundary name
std::string Face::getBoundaryName() const
{
    return boundaryName;
}

// Get face delta value
double Face::getFaceDelta() const
{
    return delta;
}

// Get face interpolation factor
double Face::getFaceInterpolationFactor() const
{
    return interpolationFactor;
}

// Get face flux
double Face::getFaceFlux() const
{
    return faceFlux;
}

// Get face zone
std::string Face::getFaceZone() const
{
    return faceZone;
}