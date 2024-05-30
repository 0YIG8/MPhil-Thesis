#ifndef DISCRETIZATION_H
#define DISCRETIZATION_H

#include "ConvectionScheme.H"
#include "Vector.H"
#include "Face.H"
#include "Cell.H"
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
void evaluateDiffusionTerm(std::vector<Face>& listOfAllFaces, ScalarField& fieldVariable, double gamma); // For constant gamma / Scalar field
void evaluateDiffusionTerm(std::vector<Face>& listOfAllFaces, ScalarField& fieldVariable, std::vector<double>& gamma); // Different gamma at every cell / Scalar field
void evaluatePositiveDiffusionTerm(std::vector<Face>& listOfAllFaces, ScalarField& fieldVariable, std::vector<double>& gamma);  // Positive Laplacian term 
void evaluateDiffusionTerm(std::vector<Face>& listOfAllFaces, VectorField& fieldVariable, double gamma); // For constant gamma / Vector field
void evaluateDiffusionTerm(std::vector<Face>& listOfAllFaces, VectorField& fieldVariable, std::vector<double>& gamma); // Different gamma at every cell / Vector field

// Convection term discretization
void evaluateConvectionTerm(std::vector<Face>& listOfAllFaces, ScalarField& fieldVariable, Vector& constVelocity, ConvectionScheme scheme); // For constant velocity
void evaluateConvectionTerm(std::vector<Face>& listOfAllFaces, VectorField& fieldVariable, Vector& constVelocity, ConvectionScheme scheme); // For constant velocity
void evaluateConvectionTerm(std::vector<Face>& listOfAllFaces, VectorField& fieldVariable, VectorField& velocityField, ConvectionScheme scheme); // For velocity field

// Flux calculation
double getFlux(Face* f, Vector constVelocity);
double getFlux(Face* f, VectorField& velocityField, int& idxP, int& idxN);
double pos(double flux);

// Under relaxation
void underRelaxImplicit(VectorField& fieldVariable, double alpha);
void underRelaxExplicit(ScalarField& fieldVariable, Eigen::VectorXd& p_old, double alpha);


#endif