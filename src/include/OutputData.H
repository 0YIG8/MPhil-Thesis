#ifndef OUTPUTDATA_H
#define OUTPUTDATA_H

#include "ScalarField.H"
#include "VectorField.H"
#include "Cell.H"
#include <vector>
#include <string>

void outputVariableDataForPlot(int nCellsx, int nCellsy, std::string filename, ScalarField fieldVariable, int iteration_number);
void outputVariableDataForPlot(int nCellsx, int nCellsy, std::string filename, VectorField fieldVariable, int iteration_number);

void outputCoeffsForPlot(VectorField fieldVariable, std::string filename, int iteration_number);
void outputCoeffsForPlot(ScalarField fieldVariable, std::string filename, int iteration_number);
void outputDataToFile(VectorField fieldVariable, std::string filename, int iteration_number, std::string matrix);
void outputDataToFile(ScalarField fieldVariable, std::string filename, int iteration_number, std::string matrix);

void outputVectorMagnitudeForPlot(int nCellsx, int nCellsy, std::string filename, VectorField fieldVariable, int iteration_number);
void outputNormalizedVectorCompForPlot(int nCellsx, int nCellsy, std::string filename, VectorField fieldVariable, int iteration_number);

void outputResidualForPlot(int nCellsx, int nCellsy, std::string filename, ScalarField fieldVariable, int iteration_number);
void outputResidualForPlot(int nCellsx, int nCellsy, std::string filename, VectorField fieldVariable, int iteration_number);
void outputAllResiduals(ScalarField p, VectorField v, std::string filename, int iteration);
void outputAllResidualNorms(ScalarField p, VectorField v, std::string filename);



#endif