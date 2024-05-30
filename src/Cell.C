#include "Cell.H"
#include "Face.H"
#include <cmath>
#include <iostream>


int Cell::nCells = 0; // Initialize no of cells to 0

/* 
    Constructor
*/

Cell::Cell(std::vector<Face *> lf)
{
    listOfFaces = lf;

    cellVolume = 0.0;
    index = nCells;
    numberOfCellNeighbours = 0;

    calcCellVolume();
    calcCellCentre();
    setFaceOwnerAndNeighbour();

    radius = cellCentre;
    omega = {0.0, 0.0, 0.0};
    cellZone = "default";

    nCells++;
}



/*
    Calculations
*/

// Calc naive cell centroid
Vector Cell::calcNaiveCellCentroid(std::vector<Face *> listOfFaces)
{
    Vector naiveCellCentroid;

    for (size_t i = 0; i < listOfFaces.size(); i++)
    {
        naiveCellCentroid = naiveCellCentroid + listOfFaces[i] -> getFaceCentre(); // Get face centroid
    }

    naiveCellCentroid = naiveCellCentroid / listOfFaces.size();

    return naiveCellCentroid;
}


// Calculate cell Centre
void Cell::calcCellCentre()
{
    Vector naiveCellCentroid = calcNaiveCellCentroid(listOfFaces); // Vector to store cell centroid

    for (size_t i = 0; i < listOfFaces.size(); i++)
    {
        Face *f = listOfFaces[i]; // Get a face
        std::vector<Point> listOfPoints = f -> getListOfPoints(); // Get its list of points
        Vector faceCentre = f->getFaceCentre(); // Get face centroid

        for (size_t j = 0; j < listOfPoints.size(); j++)
        {
            // Get 2 consecutive points & convert to vector
            Vector v1 = listOfPoints[j].pointToVector();
            Vector v2 = listOfPoints[(j + 1) % listOfPoints.size()].pointToVector();

            // Calc centroid of the tetrahedron
            Vector tetraCent = (v1 + v2 + faceCentre + naiveCellCentroid) / 4.0;

            // Calc vol of tetrahedron
            double tetraVol = fabs((v1 - faceCentre).dotMult((v2 - faceCentre).crossMult(faceCentre - naiveCellCentroid))) / 6.0;

            cellCentre = cellCentre + (tetraCent * tetraVol);
        }
    }

    // Calc cell centre
    cellCentre = cellCentre / cellVolume;
}


// Calc Cell Volume
void Cell::calcCellVolume()
{
    Vector naiveCellCentroid = calcNaiveCellCentroid(listOfFaces);

    for (size_t i = 0; i < listOfFaces.size(); i++)
    {
        Face *f = listOfFaces[i]; // Select a face
        std::vector<Point> listOfPoints = f -> getListOfPoints(); // Get list of points for this face
        Vector faceCentre = f -> getFaceCentre(); // Get centre of this face

        for (size_t j = 0; j < listOfPoints.size(); j++)
        {
            // Get 2 consecutive points & convert to vector
            Vector v1 = listOfPoints[j].pointToVector();
            Vector v2 = listOfPoints[(j + 1) % listOfPoints.size()].pointToVector();

            Vector cross = (v1 - faceCentre).crossMult(v2 - faceCentre);
            double tempVol = fabs(cross.dotMult(naiveCellCentroid - faceCentre)) / 6.0;
            cellVolume += tempVol;
        }
    }
}


// Calc radius - for MRF
void Cell::calcCellRadius(Vector origin)
{
    radius = cellCentre - origin;
}



/* 
    Set Functions
*/

// Set owner and neighbour of a cell's faces
void Cell::setFaceOwnerAndNeighbour()
{
    for (size_t i = 0; i < listOfFaces.size(); i++)
    {
        Face *f = listOfFaces[i]; // Get a face

        // Calculate dot product of vector between (face and cell centres) and face normal
        double dotProd = (f->getFaceCentre() - cellCentre).dotMult(f->getFaceNormal());
        
        // Set cell as owner or neighbour of face
        if (dotProd >= 0) // i.e. the face belongs to the owner cell
        {
            f->setFaceOwner(index);
        }
        else // If dotprod < 0 i.e. the face belongs to a neighbour cell
        {
            f->setFaceNeighbour(index);
        }
    }
}


// Set cell neighbours
void Cell::setCellNeighbours()
{
    for (size_t i = 0; i < listOfFaces.size(); i++)
    {
        if (listOfFaces[i]->getIfBoundaryFace() == true)
        {
            continue;
        }
        else if(listOfFaces[i]->getFaceOwner() == index)
        {
            listOfNeighbourCells.push_back(listOfFaces[i]->getFaceNeighbour());
            numberOfCellNeighbours++;
        }
        else
        {
            listOfNeighbourCells.push_back(listOfFaces[i]->getFaceOwner());
            numberOfCellNeighbours++;
        }
    }
}


// Set cell zone
void Cell::setCellZone(std::string zoneName)
{
    cellZone = zoneName;
}

// Set cell omega value
void Cell::setCellOmega(Vector omegaVal)
{
    omega = omegaVal;
}



/*
    Get Functions
*/

// Return Cell Centre
Vector Cell::getCellCentre() const
{
    return cellCentre;
}

// Return Cell Volume
double Cell::getCellVolume() const
{
    return cellVolume;
}

// Return Cell index
int Cell::getCellIndex() const
{
    return index;
}

// Get List of Faces
std::vector<Face *> Cell::getListOfFaces() const
{
    return listOfFaces;
}

// Get list of cell neighbours
std::vector<int> Cell::getCellNeighbours() const
{
    return listOfNeighbourCells;
}

// Get number of cell neighbours
int Cell::getNumberOfCellNeighbours() const
{
    return numberOfCellNeighbours;
}

// Get cell radius
Vector Cell::getCellRadius() const
{
    return radius;
}

// Get cell omega
Vector Cell::getCellOmega() const
{
    return omega;
}

// Get cell zone
std::string Cell::getCellZone() const
{
    return cellZone;
}