#ifndef BOUNDARYTYPE_H
#define BOUNDARYTYPE_H

#include <iostream>
#include <map>


/*
    Boundary Types
*/
enum class BoundaryType 
{
    Dirichlet,
    Neumann,
    Mixed,
    velocitySRF,
    velocityMRF,
    empty
};


// Function to convert string to BoundaryType enum
inline BoundaryType stringToBoundaryType(const std::string& BCType) 
{
    if (BCType == "Dirichlet") 
    {
        return BoundaryType::Dirichlet;
    } 
    else if (BCType == "Neumann") 
    {
        return BoundaryType::Neumann;
    } 
    else if (BCType == "Mixed") 
    {
        return BoundaryType::Mixed;
    }
    else if (BCType == "velocitySRF")
    {
        return BoundaryType::velocitySRF;
    }
    else if (BCType == "velocityMRF")
    {
        return BoundaryType::velocityMRF;
    }
    else if (BCType == "empty" || BCType == "Empty")
    {
        return BoundaryType::empty;
    }
    else 
    {
        std::cerr << "Invalid boundary type [" << BCType << "]!" << std::endl;
        exit(1);
    }
};


#endif