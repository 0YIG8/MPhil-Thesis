#include "MRFMesh.H"

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

MRFMesh::MRFMesh(Vector origin, std::string pointsFile, std::string facesFile, std::string cellsFile, std::string boundariesFile, std::string cellZonesFile, std::string faceZonesFile)
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

    // Face Zones
    if (faceZonesFile != "noInput")
    {
        // Read and Set Face Zones
        readFaceZonesFile(faceZonesFile, listOfAllFaces, listOfAllFaceZones, numberOfFacesInZone);
    }

    // Read and Set Cell Zones
    readCellZonesFile(cellZonesFile, listOfAllCells, listOfAllCellZones, numberOfCellsInZone);
}


// Default constructor
MRFMesh::MRFMesh()
{
    std::string pointsFile = "noInput";
    std::string facesFile = "noInput";
    std::string cellsFile = "noInput";
    std::string boundariesFile = "noInput";
    std::string cellZonesFile = "noInput";
    std::string faceZonesFile = "noInput";
}




/*
    Set Functions
*/
void MRFMesh::setOmegaForZone(Vector omegaVal, std::string zoneName)
{
    for (int i = 0; i < getNumberOfCellsInMesh(); i++)
    {
        if (listOfAllCells[i].getCellZone() == zoneName)
        {
            listOfAllCells[i].setCellOmega(omegaVal);
        }
    }
}



/* 
    Get Functions
*/

int MRFMesh::getNumberOfPointsInMesh()
{
    return listOfAllPoints.size();
}

int MRFMesh::getNumberOfFacesInMesh()
{
    return listOfAllFaces.size();
}

int MRFMesh::getNumberOfCellsInMesh()
{
    return listOfAllCells.size();
}

std::vector<Point> MRFMesh::getListOfPoints()
{
    return listOfAllPoints;
}

std::vector<Face> MRFMesh::getListOfFaces()
{
    return listOfAllFaces;
}

std::vector<Cell> MRFMesh::getListOfCells()
{
    return listOfAllCells;
}

int MRFMesh::getNumberOfCellZonesInMesh()
{
    return listOfAllCellZones.size();
}

std::vector<std::string> MRFMesh::getCellZoneNames()
{
    return listOfAllCellZones;
}

int MRFMesh::getNumberOfCellsInZone(std::string zoneName)
{
    auto id = std::find(listOfAllCellZones.begin(), listOfAllCellZones.end(), zoneName);

    if (id != listOfAllCellZones.end())
    {
        int index = std::distance(listOfAllCellZones.begin(), id);

        return numberOfCellsInZone[index];
    } 
    else 
    {
        std::cerr << "Cannot return number of cells in zone. Invalid cell zone name!" << std::endl;
        exit(1);
    }
}

std::vector<std::string> MRFMesh::getFaceZoneNames()
{
    return listOfAllFaceZones;
}

int MRFMesh::getNumberOfFaceZonesInMesh()
{
    return listOfAllFaceZones.size();
}

int MRFMesh::getNumberOfFacesInZone(std::string zoneName)
{
    auto id = std::find(listOfAllFaceZones.begin(), listOfAllFaceZones.end(), zoneName);

    if (id != listOfAllFaceZones.end())
    {
        int index = std::distance(listOfAllFaceZones.begin(), id);

        return numberOfFacesInZone[index];
    } 
    else 
    {
        std::cerr << "Cannot return number of faces in zone. Invalid face zone name!" << std::endl;
        exit(1);
    }
}