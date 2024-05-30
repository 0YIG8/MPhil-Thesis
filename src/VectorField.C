#include "VectorField.H"

#include <iostream>
#include <Eigen/Dense>



/*
    Constructor
*/

VectorField::VectorField(std::vector<Cell>& listOfAllCells)
{
    size = listOfAllCells.size();  // Get number of cells in mesh

    // Resize all matrices -> where, resize(nrow, ncol)
    coeffs.resize(size, size);
    sourceX.resize(size, 1);
    sourceY.resize(size, 1);
    sourceZ.resize(size, 1);
    variableX.resize(size, 1);
    variableY.resize(size, 1);
    variableZ.resize(size, 1);
    residualX.resize(size, 1);
    residualY.resize(size, 1);
    residualZ.resize(size, 1);

    // Initialize all variables to zero
    coeffs.setZero();
    sourceX.setZero();
    sourceY.setZero();
    sourceZ.setZero();
    variableX.setOnes();
    variableY.setOnes();
    variableZ.setOnes();
}



/*
    Functions
*/

// Set initial condition for the field variable -- for same value in all cells
void VectorField::setInitialCondition(double valueX, double valueY, double valueZ)
{
    variableX *= valueX;
    variableY *= valueY;
    variableZ *= valueZ;
}


// Set initial condition for field variable (-omega x radius) -- for SRF and MRF
void VectorField::setInitialCondition(std::vector<Cell>& listOfCells)
{
    if (int(listOfCells.size()) == getMatrixSize())
    {
        for (size_t i = 0; i < listOfCells.size(); i++)
        {
            // Set initial condition
            Vector vel = (listOfCells[i].getCellOmega().crossMult(listOfCells[i].getCellRadius())) * -1.0;

            variableX(i) = vel[0];
            variableY(i) = vel[1];
            variableZ(i) = vel[2];
        }
    }
    else 
    {
        std::cerr << "Mismatch between number of cells in mesh and field matrix size!" << std::endl;
        exit(1);
    }
}


// Reset coeff matrix and source matrix (eg. to be used at the end of a time step)
void VectorField::resetMatrices()
{
    coeffs.setZero();

    sourceX.setZero();
    sourceY.setZero();
    sourceZ.setZero();
}


// Set a boundary condition for the field
void VectorField::setBoundaryCondition(std::string patchName, BoundaryType BCType, double valueX, double valueY, double valueZ)
{
    hasBoundaryCondition = true;

    nameOfBoundaryPatch.push_back(patchName);
    typeOfBoundaryCondition.push_back(BCType);
    boundaryValueX.push_back(valueX);
    boundaryValueY.push_back(valueY);
    boundaryValueZ.push_back(valueZ);
}


// Solve the system
void VectorField::solveEquation()
{
    variableX = coeffs.lu().solve(sourceX);
    variableY = coeffs.lu().solve(sourceY);
    variableZ = coeffs.lu().solve(sourceZ);
}


// Add source contribution
void VectorField::addSourceContribution(int idx, double valueX, double valueY, double valueZ)
{
    sourceX(idx) += valueX;
    sourceY(idx) += valueY;
    sourceZ(idx) += valueZ;
}


// Calculate residual
void VectorField::calculateResidual()
{
    residualX = sourceX - (coeffs * variableX);
    residualY = sourceY - (coeffs * variableY);
    residualZ = sourceZ - (coeffs * variableZ);
}

// Calculate residual norm
void VectorField::calculateResidualNorm()
{
    double sumX = 0.0, sumY = 0.0, sumZ = 0.0;

    for (int i = 0; i < size; i++)
    {
        sumX += residualX(i);
        sumY += residualY(i);
        sumZ += residualZ(i);
    }

    residualNormX.push_back(sumX);
    residualNormY.push_back(sumY);
    residualNormZ.push_back(sumZ);
}


/* 
    Get Functions
*/

int VectorField::getMatrixSize()
{
    return size;
}


bool VectorField::getIfHasBoundaryCondition()
{
    return hasBoundaryCondition;
}