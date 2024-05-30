#ifndef SRFMESH_H
#define SRFMESH_H

#include <vector>
#include <string>
#include "Vector.H"
#include "Point.H"
#include "Face.H"
#include "Cell.H"

class SRFMesh
{
    private:
        std::vector<Point> listOfAllPoints;
        std::vector<Face> listOfAllFaces;
        std::vector<Cell> listOfAllCells;

        std::vector<std::vector<int>> listOfFaceIndices;
        std::vector<std::vector<int>> listOfCellIndices;
        std::vector<std::vector<int>> listOfAllBoundaries;

    public:
        SRFMesh(std::string pointsFile, std::string facesFile, std::string cellsFile, std::string boundariesFile, Vector origin);

        void setOmegaForAllCells(Vector omegaVal);

        int getNumberOfPointsInMesh();
        int getNumberOfFacesInMesh();
        int getNumberOfCellsInMesh();
        std::vector<Point> getListOfPoints();
        std::vector<Face> getListOfFaces();
        std::vector<Cell> getListOfCells();
};

#endif