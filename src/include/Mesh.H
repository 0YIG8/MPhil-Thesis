#ifndef MESH_H
#define MESH_H

#include <vector>
#include <string>
#include "Point.H"
#include "Face.H"
#include "Cell.H"

class Mesh
{
    private:
        std::vector<Point> listOfAllPoints;
        std::vector<Face> listOfAllFaces;
        std::vector<Cell> listOfAllCells;

        std::vector<std::vector<int>> listOfFaceIndices;
        std::vector<std::vector<int>> listOfCellIndices;
        std::vector<std::vector<int>> listOfAllBoundaries;

    public:
        Mesh(std::string pointsFile, std::string facesFile, std::string cellsFile, std::string boundariesFile);

        int getNumberOfPointsInMesh();
        int getNumberOfFacesInMesh();
        int getNumberOfCellsInMesh();
        std::vector<Point> getListOfPoints();
        std::vector<Face> getListOfFaces();
        std::vector<Cell> getListOfCells();
};

#endif