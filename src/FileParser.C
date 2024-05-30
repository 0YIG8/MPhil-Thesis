#include "FileParser.H"

#include "Vector.H"
#include "Point.H"
#include "Face.H"
#include "Cell.H"
#include "BoundaryPatch.H"

#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <cctype>
#include <algorithm>


/* 
    Function that reads in a points file and returns a list of points 
*/
std::vector<Point> readPointsFile(std::string fileName)
{
    std::ifstream pointsFileIn;
    pointsFileIn.open(fileName);

    int nPts; // No. of points in the file
    std::string X, Y, Z; // To get and store x, y, z in string form
    double x, y, z; // To be passed to Point in vector
    std::string lineStr; // variable to save lines in string form
    std::vector<Point> listOfPoints;

    if (pointsFileIn.is_open() && !pointsFileIn.fail()) //if the file is open
    {
        std::getline(pointsFileIn, lineStr, '\n'); // Read no. of pts
        nPts = std::stoi(lineStr); // Store no of pts
        std::getline(pointsFileIn, lineStr, '\n'); // Skip the line with the bracket

        for (int i = 0; i < nPts; i++) //while the end of file is NOT reached
        {
            std::getline(pointsFileIn, lineStr, '(');
            std::getline(pointsFileIn, X, ' ');
            std::getline(pointsFileIn, Y, ' ');
            std::getline(pointsFileIn, Z, ')');

            x = std::stod(X);
            y = std::stod(Y);
            z = std::stod(Z);

            Point p(x, y, z);
            listOfPoints.push_back(p);
        }
    }
    else // If the file is not open
    {
        std::cerr << "Unable to open " << fileName << "!" << std::endl;
        exit(1);
    }

    return listOfPoints;
}



/* 
    Function that reads a list of indices of points in a face, for all faces, and returns it
*/
std::vector<std::vector<int>> readFacesFile(std::string fileName)
{
    std::ifstream facesFileIn;
    facesFileIn.open(fileName);

    int nFcs; // No. of faces in the file
    int nPts; // No. of points for each face
    std::string lineStr; // variable to save lines in string form
    std::vector<std::vector<int>> listOfFaces;

    if (facesFileIn.is_open()) //if the file is open
    {
        std::getline(facesFileIn, lineStr, '\n'); // Read no. of faces
        nFcs = std::stoi(lineStr); // Store no of faces
        std::getline(facesFileIn, lineStr, '\n'); // Skip the line with the bracket

        for (int i = 0; i < nFcs; i++)
        {   
            // Skip whitespace and read in number of points in the face
            std::getline(facesFileIn, lineStr, '(');
            std::istringstream ss(lineStr);
            ss >> std::ws;
            std::getline(ss, lineStr, '(');
            nPts = std::stoi(lineStr);

            std::vector<int> face; // Vector to store index of points for a single face
            int pointIdx; // Index of a point in a face

            // Store index of each point in the face std::vector
            for (int j = 0; j < nPts; j++)
            {
                if (j < nPts-1)
                {
                    std::getline(facesFileIn, lineStr, ' ');
                }
                else
                {
                    std::getline(facesFileIn, lineStr, '\n');
                }
                pointIdx = std::stoi(lineStr);
                face.push_back(pointIdx);
            }

            // Add the vector to the list of faces
            listOfFaces.push_back(face);
        }
    }
    else // If file is not open
    {
        std::cerr << "Unable to open " << fileName << "!" << std::endl;
        exit(1);
    }

    return listOfFaces;
}



/* 
    Function that reads a list of indices of faces in a cell, for all cells, and returns it 
*/
std::vector<std::vector<int>> readCellsFile(std::string fileName)
{
    std::ifstream cellsFileIn;
    cellsFileIn.open(fileName);

    int nCls; // No. of cells in the file
    int nFcs; // No. of faces belonging to the cell
    std::string lineStr; // variable to save lines in string form
    std::vector<std::vector<int>> listOfCells;

    if (cellsFileIn.is_open()) //if the file is open
    {
        std::getline(cellsFileIn, lineStr, '\n'); // Read no. of cells
        nCls = std::stoi(lineStr); // Store no of cells
        std::getline(cellsFileIn, lineStr, '\n'); // Skip the line with the bracket

        for (int i = 0; i < nCls; i++)
        {   
            // Skip whitespace and read in number of faces belonging to the cell
            std::getline(cellsFileIn, lineStr, '(');
            std::istringstream ss(lineStr);
            ss >> std::ws;
            std::getline(ss, lineStr, '(');
            nFcs = std::stoi(lineStr);

            std::vector<int> cell; // Vector to store index of faces for a single cell
            int faceIdx; // Index of a face in a cell

            // Store index of each face in the cell std::vector
            for (int j = 0; j < nFcs; j++)
            {
                if (j < nFcs-1)
                {
                    std::getline(cellsFileIn, lineStr, ' ');
                }
                else
                {
                    std::getline(cellsFileIn, lineStr, '\n');
                }
                faceIdx = std::stoi(lineStr);
                cell.push_back(faceIdx);
            }

            // Add the vector to the list of faces
            listOfCells.push_back(cell);
        }
    }
    else // If file is not open
    {
        std::cerr << "Unable to open " << fileName << "!" << std::endl;
        exit(1);
    }

    return listOfCells;
}



/* 
    Function that reads a list of boundary conditions and sets them
*/
void readBoundariesFile(std::string fileName, std::vector<Face>& listOfAllFaces)
{
    std::ifstream boundariesFileIn;
    boundariesFileIn.open(fileName);
    std::string lineStr; // variable to save lines in string form
    std::string bName; // variable to save the name of the boundary

    int nBCTypes = 0; // Number of differerent boundary conditions
    int nOfBFaces = 0; // Variable to store no. of boundary faces for each type 

    if (boundariesFileIn.is_open()) //if the file is open
    {
        std::getline(boundariesFileIn, lineStr, '\n'); // Read no. of BC types
        nBCTypes = std::stoi(lineStr); // Store no of BC types
        std::getline(boundariesFileIn, lineStr, '\n'); // Skip the line with the bracket

        for (int i = 0; i < nBCTypes; i++)
        {
            // Ignore whitespace and store boundary type
            std::getline(boundariesFileIn, bName, '\n');
            bName.erase(std::remove_if(bName.begin(), bName.end(), ::isspace), bName.end());

            std::getline(boundariesFileIn, lineStr, '\n');
            lineStr.erase(std::remove_if(lineStr.begin(), lineStr.end(), ::isspace), lineStr.end());
            nOfBFaces = std::stoi(lineStr);
            std::getline(boundariesFileIn, lineStr, '\n'); // Skip the line with the bracket

            // Get list of face indices of this boundary type
            std::vector<int> lf; // List of faces

            // Skip whitespace and read in indices of faces belonging to the cell
            std::getline(boundariesFileIn, lineStr, '\n');
            auto whitespace = std::find_if(lineStr.begin(), lineStr.end(), [](unsigned char ch) {return !std::isspace(ch);});
            lineStr.erase(lineStr.begin(), whitespace);
            std::istringstream ss(lineStr);
                
            for (int j = 0; j < nOfBFaces; j++)
            {
                if (j < nOfBFaces-1)
                {
                    std::getline(ss, lineStr, ' ');
                }
                else
                {
                    std::getline(ss, lineStr, '\n');
                }

                lf.push_back(std::stoi(lineStr));
            }

            std::getline(boundariesFileIn, lineStr, '\n'); // Skip the line with the bracket

            // Set up boundaries
            BoundaryPatch b(lf, bName, listOfAllFaces);
        }
    }
    else // If file is not open
    {
        std::cerr << "Unable to open " << fileName << "!" << std::endl;
        exit(1);
    }
}



/*
    Function that reads face zones file and sets them
*/
void readFaceZonesFile(std::string fileName, std::vector<Face>& listOfAllFaces, std::vector<std::string>& listOfFaceZones, std::vector<int>& numOfFacesInZone)
{
    std::ifstream faceZonesFileIn;
    faceZonesFileIn.open(fileName);
    std::string lineStr; // variable to save lines in string form
    std::string zName;

    int nFaceZones = 0; // Number of differerent face zones
    int nFacesInZone = 0; // Variable to store no. of faces for each zone

    if (faceZonesFileIn.is_open()) //if the file is open
    {
        std::getline(faceZonesFileIn, lineStr, '\n'); // Read no. of different face zones
        nFaceZones = std::stoi(lineStr); // Store no of different face zones
        std::getline(faceZonesFileIn, lineStr, '\n'); // Skip the line with the bracket

        for (int i = 0; i < nFaceZones; i++)
        {
            // Ignore whitespace and store zone name
            std::getline(faceZonesFileIn, zName, '\n');
            zName.erase(std::remove_if(zName.begin(), zName.end(), ::isspace), zName.end());
            listOfFaceZones.push_back(zName); // Add zone name to list of zones
 
            // Skip the next 3 lines
            std::getline(faceZonesFileIn, lineStr, '\n'); // Skip the line with the bracket
            std::getline(faceZonesFileIn, lineStr, '\n'); // Skip the line with type faceZone
            std::getline(faceZonesFileIn, lineStr, '\n'); // Skip the line with faceLabels

            // Get no. of faces in this zone
            std::getline(faceZonesFileIn, lineStr, '\n'); // Read the next line with no. of faces in this zone
            lineStr.erase(std::remove_if(lineStr.begin(), lineStr.end(), ::isspace), lineStr.end());
            nFacesInZone = std::stoi(lineStr); // Store no. of faces in this zone
            numOfFacesInZone.push_back(nFacesInZone); // Add number of faces in zone to list
            std::getline(faceZonesFileIn, lineStr, '\n'); // Skip the line with the bracket

            // Read in indices of faces belonging to the cell + skip whitespaces and \n
            for (int j = 0; j < nFacesInZone; j++)
            {
                std::getline(faceZonesFileIn, lineStr, '\n');
                int faceIndex = std::stoi(lineStr); // Store face index

                // Set zone name for these faces
                for (size_t k = 0; k < listOfAllFaces.size(); k++)
                {
                    if (faceIndex == listOfAllFaces[k].getFaceIndex())
                    {
                        listOfAllFaces[k].setFaceZone(zName);
                    }
                }
            }

            // Skip the next 5 lines
            std::getline(faceZonesFileIn, lineStr, '\n'); // Skip the line with the bracket
            std::getline(faceZonesFileIn, lineStr, '\n'); // Skip the line with ;
            std::getline(faceZonesFileIn, lineStr, '\n'); // Skip the line with flipMap
            std::getline(faceZonesFileIn, lineStr, '\n'); // Skip the line with }
            std::getline(faceZonesFileIn, lineStr, '\n'); // Skip the empty line            
        }
    }
    else // If file is not open
    {
        std::cerr << "Unable to open " << fileName << "!" << std::endl;
        exit(1);
    }
}



/*
    Function that reads cell zones file and sets them
*/
void readCellZonesFile(std::string fileName, std::vector<Cell>& listOfAllCells, std::vector<std::string>& listOfCellZones, std::vector<int>& numOfCellsInZone)
{
    std::ifstream cellZonesFileIn;
    cellZonesFileIn.open(fileName);
    std::string lineStr; // variable to save lines in string form
    std::string zName;

    int nCellZones = 0; // Number of differerent cell zones
    int nCellsInZone = 0; // Variable to store no. of cells for each zone

    if (cellZonesFileIn.is_open()) //if the file is open
    {
        std::getline(cellZonesFileIn, lineStr, '\n'); // Read no. of different cell zones
        nCellZones = std::stoi(lineStr); // Store no of different cell zones
        std::getline(cellZonesFileIn, lineStr, '\n'); // Skip the line with the bracket

        for (int i = 0; i < nCellZones; i++)
        {
            // Ignore whitespace and store zone name
            std::getline(cellZonesFileIn, zName, '\n');
            zName.erase(std::remove_if(zName.begin(), zName.end(), ::isspace), zName.end());
            listOfCellZones.push_back(zName); // Add zone name to list of zones

            // Skip the next 3 lines
            std::getline(cellZonesFileIn, lineStr, '\n'); // Skip the line with the bracket
            std::getline(cellZonesFileIn, lineStr, '\n'); // Skip the line with type faceZone
            std::getline(cellZonesFileIn, lineStr, '\n'); // Skip the line with faceLabels

            // Get no. of cells in this zone
            std::getline(cellZonesFileIn, lineStr, '\n'); // Read the next line with no. of cells in this zone
            lineStr.erase(std::remove_if(lineStr.begin(), lineStr.end(), ::isspace), lineStr.end());
            nCellsInZone = std::stoi(lineStr); // Store no. of cells in this zone
            numOfCellsInZone.push_back(nCellsInZone); // Add number of faces in zone to list
            std::getline(cellZonesFileIn, lineStr, '\n'); // Skip the line with the bracket

            // Read in indices of faces belonging to the cell + skip whitespaces and \n
            for (int j = 0; j < nCellsInZone; j++)
            {
                std::getline(cellZonesFileIn, lineStr, '\n');
                int cellIndex = std::stoi(lineStr); // Store cell index

                // Set zone name for these cells
                for (size_t k = 0; k < listOfAllCells.size(); k++)
                {
                    if (cellIndex == listOfAllCells[k].getCellIndex())
                    {
                        listOfAllCells[k].setCellZone(zName);
                    }
                }
            }

            // Skip the next 5 lines
            std::getline(cellZonesFileIn, lineStr, '\n'); // Skip the line with the bracket
            std::getline(cellZonesFileIn, lineStr, '\n'); // Skip the line with ;
            std::getline(cellZonesFileIn, lineStr, '\n'); // Skip the line with }
            std::getline(cellZonesFileIn, lineStr, '\n'); // Skip the empty line            
        }
    }
    else // If file is not open
    {
        std::cerr << "Unable to open " << fileName << "!" << std::endl;
        exit(1);
    }
}