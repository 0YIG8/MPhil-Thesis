#ifndef BOUNDARYPATCH_H
#define BOUNDARYPATCH_H

#include <vector>
#include <string>
#include "Face.H"


class BoundaryPatch
{
    static int nBoundaryPatches;

    private:

    public:
        BoundaryPatch(std::vector<int> listOfBoundaryFaceIndices, std::string BCName, std::vector<Face>& listOfAllFaces);
};

#endif