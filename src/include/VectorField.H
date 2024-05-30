#ifndef VECTORFIELD_H
#define VECTORFIELD_H

#include "Cell.H"
#include "BoundaryType.H"

#include <vector>
#include <Eigen/Dense>



/*
    VectorField class
*/

class VectorField
{
    private:
        int size;
        bool hasBoundaryCondition;

    public:
        VectorField(std::vector<Cell>& listOfAllCells);

        Eigen::MatrixXd coeffs; // Matrix of dimensions dynamic x dynamic of type double
        Eigen::VectorXd variableX; // Matrix of dimensions dynamic x 1 of type double
        Eigen::VectorXd variableY;
        Eigen::VectorXd variableZ;
        Eigen::VectorXd sourceX;
        Eigen::VectorXd sourceY;
        Eigen::VectorXd sourceZ;

        Eigen::VectorXd residualX;
        Eigen::VectorXd residualY;
        Eigen::VectorXd residualZ;
        std::vector<double> residualNormX;
        std::vector<double> residualNormY;
        std::vector<double> residualNormZ;

        std::vector<std::string> nameOfBoundaryPatch;
        std::vector<BoundaryType> typeOfBoundaryCondition;
        std::vector<double> boundaryValueX;
        std::vector<double> boundaryValueY;
        std::vector<double> boundaryValueZ;

        void setInitialCondition(double valueX, double valueY, double valueZ); // For Mesh
        void setInitialCondition(std::vector<Cell>& listOfCells); // For SRF or MRF

        void resetMatrices();
        void setBoundaryCondition(std::string patchName, BoundaryType BCType, double valueX, double valueY, double valueZ);
        void solveEquation();
        void addSourceContribution(int idx, double valueX, double valueY, double valueZ);


        void calculateResidual();
        void calculateResidualNorm();


        int getMatrixSize();
        bool getIfHasBoundaryCondition();
        std::vector<double> getResidualNorm();

};

#endif