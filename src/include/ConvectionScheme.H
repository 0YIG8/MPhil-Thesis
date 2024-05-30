#ifndef CONVECTIONSCHEME_H
#define CONVECTIONSCHEME_H

#include <iostream>
#include <map>

/*
    Convection Schemes
*/
enum class ConvectionScheme 
{
    Central,
    Upwind
};


// Function to convert string to ConvectionScheme enum
inline ConvectionScheme stringToConvectionScheme(const std::string& scheme) 
{
    if (scheme == "Central" || scheme == "Central Differencing" || scheme == "CD") 
    {
        return ConvectionScheme::Central;
    } 
    else if (scheme == "Upwind" || scheme == "Upwind Differencing" || scheme == "UD") 
    {
        return ConvectionScheme::Upwind;
    }
    else 
    {
        std::cerr << "Invalid convection scheme [" << scheme << "]!" << std::endl;
        exit(1);
    }
};


#endif