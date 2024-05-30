#ifndef FILEPARSER_H
#define FILEPARSER_H

#include <vector>
#include <string>
#include "Point.H"
#include "Face.H"
#include "Cell.H"

std::vector<Point> readPointsFile(std::string fileName);
std::vector<std::vector<int>> readFacesFile(std::string fileName);
std::vector<std::vector<int>> readCellsFile(std::string fileName);
void readBoundariesFile(std::string fileName, std::vector<Face>& listOfAllFaces);

void readFaceZonesFile(std::string fileName, std::vector<Face>& listOfAllFaces, std::vector<std::string>& listOfFaceZones, std::vector<int>& numOfFacesInZone);
void readCellZonesFile(std::string fileName, std::vector<Cell>& listOfAllCells, std::vector<std::string>& listOfCellZones, std::vector<int>& numOfCellsInZone);


#endif