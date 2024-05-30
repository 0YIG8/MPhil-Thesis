#ifndef MRFMESH_H
#define MRFMESH_H

#include <vector>
#include <string>
#include "Vector.H"
#include "Point.H"
#include "Face.H"
#include "Cell.H"

class MRFMesh
{
    private:
        std::vector<Point> listOfAllPoints;
        std::vector<Face> listOfAllFaces;
        std::vector<Cell> listOfAllCells;

        std::vector<std::string> listOfAllFaceZones;
        std::vector<int> numberOfFacesInZone;

        std::vector<std::string> listOfAllCellZones;
        std::vector<int> numberOfCellsInZone;

        std::vector<std::vector<int>> listOfFaceIndices;
        std::vector<std::vector<int>> listOfCellIndices;
        std::vector<std::vector<int>> listOfAllBoundaries;

    public:
        MRFMesh(Vector origin, std::string pointsFile, std::string facesFile, std::string cellsFile, std::string boundariesFile, std::string cellZonesFile, std::string faceZonesFile = "noInput");
        MRFMesh();  // Default constructor

        void setOmegaForZone(Vector omegaVal, std::string zoneName);

        int getNumberOfPointsInMesh();
        int getNumberOfFacesInMesh();
        int getNumberOfCellsInMesh();
        std::vector<Point> getListOfPoints();
        std::vector<Face> getListOfFaces();
        std::vector<Cell> getListOfCells();

        std::vector<std::string> getCellZoneNames();
        int getNumberOfCellZonesInMesh();
        int getNumberOfCellsInZone(std::string zoneName);
        std::vector<std::string> getFaceZoneNames();
        int getNumberOfFaceZonesInMesh();
        int getNumberOfFacesInZone(std::string zoneName);
};

#endif