#ifndef SIMPLE_H
#define SIMPLE_H

#include "MRFMesh.H"
#include "SRFMesh.H"
#include "VectorField.H"
#include "ScalarField.H"
#include <libconfig.h++>

using namespace libconfig;

void evaluateSIMPLE(MRFMesh& mesh, VectorField& velocity, ScalarField& pressure, const Setting& root);
void evaluateSIMPLE(SRFMesh& mesh, VectorField& velocity, ScalarField& pressure, const Setting& root);


#endif