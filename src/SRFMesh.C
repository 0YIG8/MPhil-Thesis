#include "SRFMesh.H"

#include "Vector.H"
#include "Point.H"
#include "Face.H"
#include "Cell.H"
#include "FileParser.H"
#include "Discretization.H"
#include <iostream>


/* 
    Constructor
*/

SRFMesh::SRFMesh(std::string pointsFile, std::string facesFile, std::string cellsFile, std::string boundariesFile, Vector origin)
{
    // Get and store list of points from points file
    listOfAllPoints = readPointsFile(pointsFile);

    // Get and store list of faces from faces file
    listOfFaceIndices = readFacesFile(facesFile);
    for (size_t i = 0; i < listOfFaceIndices.size(); i++)
    {
        std::vector<int> m = listOfFaceIndices[i];
        std::vector<Point> listOfPoints;
        for (size_t idx = 0; idx < listOfFaceIndices[i].size(); idx++)
        {
            Point pt = Point::getPointFromIndex(m[idx], listOfAllPoints);
            listOfPoints.push_back(pt);
        }

        Face f(listOfPoints);
        listOfAllFaces.push_back(f);
    }

    // Get and store list of cells from cells file
    listOfCellIndices = readCellsFile(cellsFile);
    for (size_t i = 0; i < listOfCellIndices.size(); i++)
    {
        std::vector<int> n = listOfCellIndices[i];
        std::vector<Face *> listOfFaces;
        for (size_t idx = 0; idx < listOfCellIndices[i].size(); idx++)
        {
            int k = n[idx];
            listOfFaces.push_back(&listOfAllFaces[k]);
        }

        Cell c(listOfFaces);
        listOfAllCells.push_back(c);
    }

    // Read and Set Boundary Conditions
    readBoundariesFile(boundariesFile, listOfAllFaces);

    // Set cell neighbours, calc radius
    for (size_t i = 0; i < listOfAllCells.size(); i++)
    {
        listOfAllCells[i].setCellNeighbours();
        listOfAllCells[i].calcCellRadius(origin);
    }

    // Calc and store delta and interpolation factor values for each face
    for (size_t i = 0; i < listOfAllFaces.size(); i++)
    {
        Face* f = &listOfAllFaces[i];
        f->setFaceDelta(delta(f, listOfAllCells));
        f->setFaceInterpolationFactor(interpolationFactor(f, listOfAllCells));
    }
}



/*
    Set Functions
*/
// Set omega values for all cells
void SRFMesh::setOmegaForAllCells(Vector omegaVal)
{
    for (int i = 0; i < getNumberOfCellsInMesh(); i++)
    {
        listOfAllCells[i].setCellOmega(omegaVal);
    }
}



/* 
    Get Functions
*/

int SRFMesh::getNumberOfPointsInMesh()
{
    return listOfAllPoints.size();
}

int SRFMesh::getNumberOfFacesInMesh()
{
    return listOfAllFaces.size();
}

int SRFMesh::getNumberOfCellsInMesh()
{
    return listOfAllCells.size();
}

std::vector<Point> SRFMesh::getListOfPoints()
{
    return listOfAllPoints;
}

std::vector<Face> SRFMesh::getListOfFaces()
{
    return listOfAllFaces;
}

std::vector<Cell> SRFMesh::getListOfCells()
{
    return listOfAllCells;
}