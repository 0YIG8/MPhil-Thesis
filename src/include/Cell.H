#ifndef CELL_H
#define CELL_H

#include "Face.H"

class Cell
{
    static int nCells;  // Number of Cells

    private:
        std::vector<Face *> listOfFaces;

        double cellVolume;
        Vector cellCentre;
        int index;
        std::vector<int> listOfNeighbourCells;
        int numberOfCellNeighbours;

        Vector radius;
        Vector omega;
        std::string cellZone;

    public:
        Cell(std::vector<Face *> listOfFaces);

        Vector calcNaiveCellCentroid(std::vector<Face *> listOfFaces);
        void calcCellVolume();
        void calcCellCentre();
        void calcCellRadius(Vector origin);  

        void setFaceOwnerAndNeighbour(); // For faces of cells
        void setCellNeighbours(); 

        void setCellOmega(Vector omegaVal);
        void setCellZone(std::string zoneName);

        double getCellVolume() const;
        Vector getCellCentre() const;
        int getCellIndex() const;
        std::vector<int> getCellNeighbours() const;
        std::vector<Face *> getListOfFaces() const;
        int getNumberOfCellNeighbours() const;

        Vector getCellRadius() const;
        Vector getCellOmega() const;
        std::string getCellZone() const;
};

#endif