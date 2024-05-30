/* 
    MRF PROJECT
*/

#include "BoundaryType.H"
#include "Vector.H"
#include "Cell.H"
#include "MRFMesh.H"
#include "MovingMeshDiscretization.H"
#include "SIMPLE.H"
#include "OutputVTK.H"

#include <chrono>
#include <cstddef>
#include <string>
#include <iostream>
#include <iomanip>
#include <libconfig.h++>

using namespace libconfig;


int main(int argc, char *argv[])
{
    /*
        RECORD SIMULATION START TIME
    */
    auto tStart = std::chrono::steady_clock::now();


    /*
        SIMULATION SETUP - FROM CONFIG FILE DATA
    */
    Config cfg;
    cfg.readFile(argv[1]);
    const Setting& root = cfg.getRoot();

    const Setting& simulationSetup = root["simulationSetup"];
    // const Setting& inputParams = root["inputParams"];
    const Setting& fileLocations = root["fileLocations"];
    const Setting& velocityBCs = root["velocityBCs"];
    const Setting& pressureBCs = root["pressureBCs"];

    bool isMRF = simulationSetup.lookup("isMRF");

    // Get all vectors from config file
    Vector origin, omega;
    std::vector<double> originVect, omegaVect, initialVelocity;
    const Setting& orginSett = simulationSetup["origin"];
    const Setting& omegaSett = simulationSetup["omega"];
    const Setting& initVelSett = simulationSetup["initialVelocity"];

    for (int i = 0; i < 3; i++) 
    {
        originVect.push_back(orginSett[i]);
        omegaVect.push_back(omegaSett[i]);
        initialVelocity.push_back(initVelSett[i]);
    }
    origin = {originVect[0], originVect[1], originVect[2]};
    omega = {omegaVect[0], omegaVect[1], omegaVect[2]};


    // Set up moving mesh
    MRFMesh MRF_mesh;
    
    if (fileLocations.exists("faceZones")) // If faceZones file is provided
    {
        MRF_mesh = MRFMesh(origin, fileLocations.lookup("points"), fileLocations.lookup("faces"), fileLocations.lookup("cells"), fileLocations.lookup("boundaries"), fileLocations.lookup("cellZones"), fileLocations.lookup("faceZones"));
    }
    else // if faceZones file is not provided
    {
        MRF_mesh = MRFMesh(origin, fileLocations.lookup("points"), fileLocations.lookup("faces"), fileLocations.lookup("cells"), fileLocations.lookup("boundaries"), fileLocations.lookup("cellZones"));
    }
    

    // Set omega for moving zone cells
    MRF_mesh.setOmegaForZone(omega, "rotor");

    // Get list of cells, faces, points
    std::vector<Cell> myListOfAllCls = MRF_mesh.getListOfCells(); // Get list of cells in mesh
    std::vector<Face> myListOfAllFcs = MRF_mesh.getListOfFaces(); // Get list of all faces in the mesh
    // std::vector<Point> myListOfAllPts = MRF_mesh.getListOfPoints(); // Get list of all points in the mesh

    // Setup field variables
    // 1. velocity field
    VectorField velocity(myListOfAllCls);

    // Set initial condition for velocity
    velocity.setInitialCondition(myListOfAllCls);

    // Set velocity BCs
    for (int i = 0; i < velocityBCs.getLength(); i++)
    {
        // Get velocity BC values
        const Setting& vBC = velocityBCs[i];

        std::string patchName = vBC["patchName"];
        BoundaryType BCType = stringToBoundaryType(vBC["BCType"]);

        std::vector<double> value;
        const Setting& valueSett = vBC["value"];
        for (int j = 0; j < 3; j++) 
        {
            value.push_back(valueSett[j]);
        }

        // Set boundary condition
        velocity.setBoundaryCondition(patchName, BCType, value[0], value[1], value[2]);
    }

    // Initialize flux
    initializeFlux(myListOfAllFcs, myListOfAllCls, velocity, omega, origin, isMRF);

    // 2. Pressure field
    ScalarField pressure(myListOfAllCls);  

    // Set initial condition for pressure
    pressure.setInitialCondition(simulationSetup.lookup("initialPressure"));

    // Set pressure BCs
    for (int i = 0; i < pressureBCs.getLength(); i++)
    {
        // Get pressure BC values
        const Setting& pBC = pressureBCs[i];

        std::string patchName = pBC["patchName"];
        BoundaryType BCType = stringToBoundaryType(pBC["BCType"]);
        double value = pBC["value"];

        // Set boundary condition
        pressure.setBoundaryCondition(patchName, BCType, value);
    }



    /*
        MAIN CALC - Eg. SIMPLE
    */
    evaluateSIMPLE(MRF_mesh, velocity, pressure, root);
    std::cout << "Done!" << std::endl;



    /*
        OUTPUT SIMULATION TIME
    */
    auto tEnd = std::chrono::steady_clock::now();
    std::chrono::duration<double> duration {tEnd - tStart};
    double simTime = duration.count();

    if (simTime <= 59.9) // If simulation time was less than a minute
    {
        std::cout << "Simulation complete! Time taken = " << std::setprecision(5) << simTime << " secs." << std::endl;
    }
    else // If simulation time is more than a minute
    {
        std::cout << "Simulation complete! Time taken = " << std::setprecision(5) << simTime / 60.0 << " mins." << std::endl;
    }
   
    return 0;
}