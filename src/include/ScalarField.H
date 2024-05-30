#ifndef SCALARFIELD_H
#define SCALARFIELD_H

#include "Cell.H"
#include "BoundaryType.H"

#include <vector>
#include <string>
#include <Eigen/Dense>

class ScalarField
{
    private:
        int size;
        bool hasBoundaryCondition;

    public:
        ScalarField(std::vector<Cell>& listOfAllCells);

        Eigen::MatrixXd coeffs; // Matrix of dimensions dynamic x dynamic of type double
        Eigen::VectorXd variable; // Matrix of dimensions dynamic x 1 of type double
        Eigen::VectorXd source;

        Eigen::VectorXd residual;

        std::vector<std::string> nameOfBoundaryPatch;
        std::vector<BoundaryType> typeOfBoundaryCondition;
        std::vector<double> boundaryValue;

        std::vector<double> residualNorm;

        void setInitialCondition(double value);
        void resetMatrices();
        void setBoundaryCondition(std::string patchName, BoundaryType BCType, double value);
        void solveEquation();
        void addSourceContribution(int idx, double value);

        int getMatrixSize();
        bool getIfHasBoundaryCondition();
        double getResidualNorm();

        void calculateResidual();
        void calculateResidualNorm();

        void multiplyBy(double value, std::string matrix);
};

#endif