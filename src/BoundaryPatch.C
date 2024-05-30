#include "BoundaryPatch.H"
#include "Face.H"
#include <iostream>
#include <string>


int BoundaryPatch::nBoundaryPatches = 0;

/* 
    Constructor
*/

BoundaryPatch::BoundaryPatch(std::vector<int> listOfBoundaryFaceIndices, std::string BCName, std::vector<Face>& listOfAllFaces)
{
    // Set as boundary faces
    for (size_t i = 0; i < listOfBoundaryFaceIndices.size(); i++)
    {
        for (size_t j = 0; j < listOfAllFaces.size(); j++)
        {
            if (listOfAllFaces[j].getFaceIndex() == listOfBoundaryFaceIndices[i])
            {
                listOfAllFaces[j].setAsBoundaryFace();
                listOfAllFaces[j].setBoundaryName(BCName);
            }
        }
    }

    nBoundaryPatches++;
}