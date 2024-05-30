#ifndef MOVINGMESHDISCRETIZATION_H
#define MOVINGMESHDISCRETIZATION_H

#include "ConvectionScheme.H"
#include "Vector.H"
#include "Face.H"
#include "ScalarField.H"
#include "VectorField.H"

#include <vector>


// Discretization essentials 
double delta(Face* f, std::vector<Cell>& listOfAllCells);
double interpolationFactor(Face* f, std::vector<Cell>& listOfAllCells);

// Gradient calculation
Vector getGradientGauss(ScalarField& fieldVariable, int cellIdx, std::vector<Cell>& listOfAllCells);

// Temporal derivative discretization
void evaluateTimeDerivative(std::vector<Cell>& listOfAllCells, ScalarField& fieldVariable, double delta_t);
void evaluateTimeDerivative(std::vector<Cell>& listOfAllCells, VectorField& fieldVariable, double delta_t);

// Diffusion term (Laplacian) discretization
void evaluatePositiveDiffusionTerm(std::vector<Face>& listOfAllFaces, ScalarField& fieldVariable, std::vector<double>& gamma);  // Positive Laplacian term 
void evaluateDiffusionTerm(std::vector<Face>& listOfAllFaces, VectorField& fieldVariable, double gamma, Vector omega, Vector origin);

// Convection term discretization
void evaluateConvectionTerm(std::vector<Face>& listOfAllFaces, VectorField& fieldVariable, VectorField& velocityField, ConvectionScheme scheme, Vector omega, Vector origin); // For velocity field

// Flux calculation
double getFlux(Face* f, std::vector<Cell>& listOfAllCells, VectorField& velocityField, int& idxP, int& idxN, Vector omega, Vector origin, bool isMRF);
double pos(double flux);
void initializeFlux(std::vector<Face>& listOfFaces, std::vector<Cell>& listOfAllCells, VectorField velocity, Vector omega, Vector origin, bool isMRF);

// Under relaxation
void underRelaxImplicit(VectorField& fieldVariable, double alpha);
void underRelaxExplicit(ScalarField& fieldVariable, Eigen::VectorXd& p_old, double alpha);

#endif