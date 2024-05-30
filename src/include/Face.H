#ifndef FACE_H
#define FACE_H

#include <vector>
#include "Point.H"


class Face
{
    static int nFaces;  // Number of faces

    private:
        std::vector<Point> listOfPoints;

        Vector faceCentre;
        Vector faceNormal;
        double faceArea;
        int index;
        int ownerIndex;
        int neighbourIndex;
        double delta;
        double interpolationFactor;
        double faceFlux;

        bool isBoundaryFace;
        std::string boundaryName;

        bool isMovingFace;
        std::string faceZone;

    public:
        Face(std::vector<Point> listOfPoints);

        void calcFaceArea();
        void calcFaceNormal();
        void calcFaceCentre();

        void setFaceOwner(int index); // When the face belongs to owner cell
        void setFaceNeighbour(int index); // When the face belongs to neighbour cell
        void setAsBoundaryFace();
        void setBoundaryName(std::string boundaryName);
        void setFaceDelta(double deltaValue);
        void setFaceInterpolationFactor(double intFactValue);
        void setFaceFlux(double flux);
        void setAsMovingFace(); // For MRF
        void setFaceZone(std::string zoneName);

        double getFaceArea() const;
        Vector getFaceNormal() const;
        Vector getFaceCentre() const;
        int getFaceIndex() const;
        int getFaceOwner() const;
        int getFaceNeighbour() const;
        std::vector<Point> getListOfPoints() const;
        bool getIfBoundaryFace() const;
        std::string getBoundaryName() const;
        double getFaceDelta() const;
        double getFaceInterpolationFactor() const;
        double getFaceFlux() const;
        std::string getFaceZone() const;

};

#endif