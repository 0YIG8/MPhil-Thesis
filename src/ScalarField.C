#include "ScalarField.H"
#include <Eigen/Dense>
#include <iostream>
#include <fstream>
#include <sstream>


/*
    Constructor
*/

ScalarField::ScalarField(std::vector<Cell>& listOfAllCells)
{
    size = listOfAllCells.size();  // Get number of cells in mesh
    hasBoundaryCondition = false;

    // Resize all matrices -> where, resize(nrow, ncol)
    coeffs.resize(size, size);
    source.resize(size, 1);
    variable.resize(size, 1);
    residual.resize(size, 1);

    // Initialize all variables to zero
    coeffs.setZero();
    source.setZero();
    variable.setOnes();
}



/*
    Functions
*/

// Set initial condition for the field variable -- for same value in all cells
void ScalarField::setInitialCondition(double value)
{
    variable *= value;
}


// Reset coeff matrix and source matrix (eg. to be used at the end of a time step)
void ScalarField::resetMatrices()
{
    coeffs.setZero();
    source.setZero();
}


// Set a boundary condition for the field
void ScalarField::setBoundaryCondition(std::string patchName, BoundaryType BCType, double value)
{
    hasBoundaryCondition = true;

    nameOfBoundaryPatch.push_back(patchName);
    typeOfBoundaryCondition.push_back(BCType);
    boundaryValue.push_back(value);
}


// Solve the system
void ScalarField::solveEquation()
{
    variable = coeffs.fullPivHouseholderQr().solve(source);
}


// Add source contribution
void ScalarField::addSourceContribution(int idx, double value)
{
    source(idx) += value;
}


// Calculate residual
void ScalarField::calculateResidual()
{
    // for (int i = 0; i < size; i++)
    // {
    residual = source - (coeffs * variable);
    // }
}

// Calculate residual norm
void ScalarField::calculateResidualNorm()
{
    double sum = 0.0;

    for (int i = 0; i < size; i++)
    {
        sum += residual(i);
    }

    residualNorm.push_back(sum);
}



/* 
    Get Functions
*/

int ScalarField::getMatrixSize()
{
    return size;
}


bool ScalarField::getIfHasBoundaryCondition()
{
    return hasBoundaryCondition;
}


// For negative terms
void ScalarField::multiplyBy(double value, std::string matrix)
{
    if (matrix == "coeffs")
    {
        coeffs * value;
    }
    else if (matrix == "source")
    {
        source * value;
    }
    else if (matrix == "variable")
    {
        variable * value;
    }
}