#include "Mesh.H"
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

Mesh::Mesh(std::string pointsFile, std::string facesFile, std::string cellsFile, std::string boundariesFile)
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

    // Set cell neighbours
    for (size_t i = 0; i < listOfAllCells.size(); i++)
    {
        listOfAllCells[i].setCellNeighbours();
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
    Get Functions
*/

int Mesh::getNumberOfPointsInMesh()
{
    return listOfAllPoints.size();
}

int Mesh::getNumberOfFacesInMesh()
{
    return listOfAllFaces.size();
}

int Mesh::getNumberOfCellsInMesh()
{
    return listOfAllCells.size();
}

std::vector<Point> Mesh::getListOfPoints()
{
    return listOfAllPoints;
}

std::vector<Face> Mesh::getListOfFaces()
{
    return listOfAllFaces;
}

std::vector<Cell> Mesh::getListOfCells()
{
    return listOfAllCells;
}