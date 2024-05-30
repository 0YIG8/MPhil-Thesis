#ifndef OUTPUTVTK_H
#define OUTPUTVTK_H

#include "Point.H"
#include "Cell.H"
#include "ScalarField.H"
#include "VectorField.H"
#include "MRFMesh.H"
#include "SRFMesh.H"

void outputVTK_CavityData(std::vector<Point>& listOfAllPoints, ScalarField& pressure, VectorField& velocity, std::vector<int> nCells, int iteration, std::string filename);
void outputVTK_MRFData(MRFMesh& mesh, ScalarField& pressure, VectorField& velocity, int iteration, std::string filename);
void outputVTK_SRFData(SRFMesh& mesh, ScalarField& pressure, VectorField& velocity, int iteration, std::string filename);

#endif