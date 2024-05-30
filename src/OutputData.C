#include "OutputData.H"

#include <iomanip>
#include <iostream> 
#include <fstream>
#include "ScalarField.H"
#include "VectorField.H"
#include <string>
#include <Eigen/Dense>


/*
    FUNCTIONS TO OUTPUT DATA
*/

// Output variable data for Scalar field
void outputVariableDataForPlot(int nCellsx, int nCellsy, std::string filename, ScalarField fieldVariable, int iteration_number)
{
    // Output final time step data to a single file
    std::ofstream output(filename + "_" + std::to_string(iteration_number) + ".txt");

    output << "y_coord";
    for (int i = 0; i < nCellsx; i++)
    {
        output << "," << i;
    }
    output << std::endl;

    for (int z = 0; z < nCellsy; z++)
    {
        output << z;
        for (int i = z*nCellsx; i < (z+1)*nCellsx; i++)
        {
            output << "," << fieldVariable.variable(i);
        }
        output << std::endl;
    }
}


// Output variable data for Vector field
void outputVariableDataForPlot(int nCellsx, int nCellsy, std::string filename, VectorField fieldVariable, int iteration_number)
{
    // x output
    std::ofstream outputX (filename + "_X_" + std::to_string(iteration_number) + ".txt");

    outputX << "y_coord";
    for (int i = 0; i < nCellsx; i++)
    {
        outputX << "," << i;
    }
    outputX << std::endl;

    for (int z = 0; z < nCellsy; z++)
    {
        outputX << z;
        for (int i = z*nCellsx; i < (z+1)*nCellsx; i++)
        {
            outputX << "," << fieldVariable.variableX(i);
        }
        outputX << std::endl;
    }


    // y output
    std::ofstream outputY (filename + "_Y_" + std::to_string(iteration_number) + ".txt");

    outputY << "y_coord";
    for (int i = 0; i < nCellsx; i++)
    {
        outputY << "," << i;
    }
    outputY << std::endl;

    for (int z = 0; z < nCellsy; z++)
    {
        outputY << z;
        for (int i = z*nCellsx; i < (z+1)*nCellsx; i++)
        {
            outputY << "," << fieldVariable.variableY(i);
        }
        outputY << std::endl;
    }


    // z output
    std::ofstream outputZ (filename + "_Z_" + std::to_string(iteration_number) + ".txt");

    outputZ << "y_coord";
    for (int i = 0; i < nCellsx; i++)
    {
        outputZ << "," << i;
    }
    outputZ << std::endl;

    for (int z = 0; z < nCellsy; z++)
    {
        outputZ << z;
        for (int i = z*nCellsx; i < (z+1)*nCellsx; i++)
        {
            outputZ << "," << fieldVariable.variableZ(i);
        }
        outputZ << std::endl;
    }
}


// Output variable magnitude for a vector field
void outputVectorMagnitudeForPlot(int nCellsx, int nCellsy, std::string filename, VectorField fieldVariable, int iteration_number)
{
    std::ofstream outputX (filename + "_Magnitude_" + std::to_string(iteration_number) + ".txt");

    outputX << "y_coord";
    for (int i = 0; i < nCellsx; i++)
    {
        outputX << "," << i;
    }
    outputX << std::endl;

    for (int z = 0; z < nCellsy; z++)
    {
        outputX << z;
        for (int i = z*nCellsx; i < (z+1)*nCellsx; i++)
        {
            outputX << "," << sqrt(fieldVariable.variableX(i) * fieldVariable.variableX(i) + fieldVariable.variableY(i) * fieldVariable.variableY(i) + fieldVariable.variableZ(i) * fieldVariable.variableZ(i));
        }
        outputX << std::endl;
    }
}


// Output normalised vector components for plot
void outputNormalizedVectorCompForPlot(int nCellsx, int nCellsy, std::string filename, VectorField fieldVariable, int iteration_number)
{
    // x output
    std::ofstream outputX (filename + "_X_Normalized_" + std::to_string(iteration_number) + ".txt");

    outputX << "y_coord";
    for (int i = 0; i < nCellsx; i++)
    {
        outputX << "," << i;
    }
    outputX << std::endl;

    for (int z = 0; z < nCellsy; z++)
    {
        outputX << z;
        for (int i = z*nCellsx; i < (z+1)*nCellsx; i++)
        {
            outputX << "," << fieldVariable.variableX(i) / std::sqrt(fieldVariable.variableX(i) * fieldVariable.variableX(i) + fieldVariable.variableY(i) * fieldVariable.variableY(i));
        }
        outputX << std::endl;
    }


    // y output
    std::ofstream outputY (filename + "_Y_Normalized_" + std::to_string(iteration_number) + ".txt");

    outputY << "y_coord";
    for (int i = 0; i < nCellsx; i++)
    {
        outputY << "," << i;
    }
    outputY << std::endl;

    for (int z = 0; z < nCellsy; z++)
    {
        outputY << z;
        for (int i = z*nCellsx; i < (z+1)*nCellsx; i++)
        {
            outputY << "," << fieldVariable.variableY(i) / std::sqrt(fieldVariable.variableX(i) * fieldVariable.variableX(i) + fieldVariable.variableY(i) * fieldVariable.variableY(i));
        }
        outputY << std::endl;
    }
}


// Output coeff matrix for a vector field
void outputCoeffsForPlot(VectorField fieldVariable, std::string filename, int iteration_number)
{
    std::ofstream output(filename + "Coeffs_" + std::to_string(iteration_number) + ".txt");

    int size = fieldVariable.getMatrixSize();

    output << "row-col";
    for (int i = 0; i < size; i++)
    {
        output << "," << i;
    }
    output << std::endl;

    for (int i = 0; i < size; i++)
    {
        output << i;
        for (int j = 0; j < size; j++)
        {
            output << "," << fieldVariable.coeffs(i, j);
        }
        output << std::endl;
    }
}


// Output coeff data for a scalar field
void outputCoeffsForPlot(ScalarField fieldVariable, std::string filename, int iteration_number)
{
    // Output final time step data to a single file
    std::ofstream output(filename + "Coeffs_" + std::to_string(iteration_number) + ".txt");

    int size = fieldVariable.getMatrixSize();

    output << "row-col";
    for (int i = 0; i < size; i++)
    {
        output << "," << i;
    }
    output << std::endl;

    for (int i = 0; i < size; i++)
    {
        output << i;
        for (int j = 0; j < size; j++)
        {
            output << "," << fieldVariable.coeffs(i, j);
        }
        output << std::endl;
    }
}


// Output source term data for a vector field
void outputDataToFile(VectorField fieldVariable, std::string filename, int iteration_number, std::string matrix)
{
    if (matrix == "Source" || matrix == "source")
    {
        std::ofstream output(filename + "Source_" + std::to_string(iteration_number) + ".txt");

        int size = fieldVariable.getMatrixSize();

        output << "idx  X   Y   Z" << std::endl;

        for (int i = 0; i < size; i++)
        {
            output << i << "    " << fieldVariable.sourceX(i) << "  " << fieldVariable.sourceY(i) << "  " << fieldVariable.sourceZ(i) << std::endl;
        }
    }
    else if (matrix == "Variable" || matrix == "variable")
    {
        std::ofstream output(filename + "Variable_" + std::to_string(iteration_number) + ".txt");

        int size = fieldVariable.getMatrixSize();

        output << "idx  X   Y   Z" << std::endl;

        for (int i = 0; i < size; i++)
        {
            output << i << "    " << fieldVariable.variableX(i) << "  " << fieldVariable.variableY(i) << "  " << fieldVariable.variableZ(i) << std::endl;
        }
    }
}


// Output source term data for a scalar field
void outputDataToFile(ScalarField fieldVariable, std::string filename, int iteration_number, std::string matrix)
{
    if (matrix == "Source" || matrix == "source")
    {
        std::ofstream output(filename + "Source_" + std::to_string(iteration_number) + ".txt");

        int size = fieldVariable.getMatrixSize();

        output << "idx  X" << std::endl;

        for (int i = 0; i < size; i++)
        {
            output << i << "    " << fieldVariable.source(i) << std::endl;
        }
    }
    else if (matrix == "Variable" || matrix == "variable")
    {
        std::ofstream output(filename + "Variable_" + std::to_string(iteration_number) + ".txt");

        int size = fieldVariable.getMatrixSize();

        output << "idx  X" << std::endl;

        for (int i = 0; i < size; i++)
        {
            output << i << "    " << fieldVariable.variable(i) << std::endl;
        }
    }
}



// Output residual data for Scalar field
void outputResidualForPlot(int nCellsx, int nCellsy, std::string filename, ScalarField fieldVariable, int iteration_number)
{
    // Output final time step data to a single file
    std::ofstream output(filename + "_Residual_" + std::to_string(iteration_number) + ".txt");

    output << "y_coord";
    for (int i = 0; i < nCellsx; i++)
    {
        output << "," << i;
    }
    output << std::endl;

    for (int z = 0; z < nCellsy; z++)
    {
        output << z;
        for (int i = z*nCellsx; i < (z+1)*nCellsx; i++)
        {
            output << std::setprecision(100) << "," << fieldVariable.residual(i);
        }
        output << std::endl;
    }
}


// Output residuals for Vector field
void outputResidualForPlot(int nCellsx, int nCellsy, std::string filename, VectorField fieldVariable, int iteration_number)
{
    // x output
    std::ofstream outputX (filename + "_ResidualX_" + std::to_string(iteration_number) + ".txt");

    outputX << "y_coord";
    for (int i = 0; i < nCellsx; i++)
    {
        outputX << "," << i;
    }
    outputX << std::endl;

    for (int z = 0; z < nCellsy; z++)
    {
        outputX << z;
        for (int i = z*nCellsx; i < (z+1)*nCellsx; i++)
        {
            outputX << std::setprecision(100) << "," << fieldVariable.residualX(i);
        }
        outputX << std::endl;
    }


    // y output
    std::ofstream outputY (filename + "_ResidualY_" + std::to_string(iteration_number) + ".txt");

    outputY << "y_coord";
    for (int i = 0; i < nCellsx; i++)
    {
        outputY << "," << i;
    }
    outputY << std::endl;

    for (int z = 0; z < nCellsy; z++)
    {
        outputY << z;
        for (int i = z*nCellsx; i < (z+1)*nCellsx; i++)
        {
            outputY << std::setprecision(100) << "," << fieldVariable.residualY(i);
        }
        outputY << std::endl;
    }


    // z output
    std::ofstream outputZ (filename + "_ResidualZ_" + std::to_string(iteration_number) + ".txt");

    outputZ << "y_coord";
    for (int i = 0; i < nCellsx; i++)
    {
        outputZ << "," << i;
    }
    outputZ << std::endl;

    for (int z = 0; z < nCellsy; z++)
    {
        outputZ << z;
        for (int i = z*nCellsx; i < (z+1)*nCellsx; i++)
        {
            outputZ << std::setprecision(100) << "," << fieldVariable.residualZ(i);
        }
        outputZ << std::endl;
    }
}


void outputAllResiduals(ScalarField p, VectorField v, std::string filename, int iteration)
{
    int size = p.getMatrixSize();

    std::ofstream output(filename + "_ResLinePlot_" + std::to_string(iteration) + ".txt");

    output << "Cell_index,P_res,Vx_res,Vy_res,Vz_res" << std::endl;

    for (int z = 0; z < size; z++)
    {
        output << std::setprecision(100) << z << "," << p.residual(z) << "," << v.residualX(z) << "," << v.residualY(z) << "," << v.residualZ(z) << std::endl;
    }
}

// To be used after (for eg) the SIMPLE loop
void outputAllResidualNorms(ScalarField p, VectorField v, std::string filename)
{
    std::ofstream output(filename + "ResidualNorm.txt");

    output << "Iteration,P_resNorm,Vx_resNorm,Vy_resNorm,Vz_resNorm" << std::endl;

    for (size_t z = 0; z < p.residualNorm.size(); z++)
    {
        output << std::setprecision(10) << z << "," << p.residualNorm[z] << "," << v.residualNormX[z] << "," << v.residualNormY[z] << "," << v.residualNormZ[z] << std::endl;
    }
}