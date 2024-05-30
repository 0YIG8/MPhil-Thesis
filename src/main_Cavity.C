/* 
    LID-DRIVEN CAVITY WITH VTK OUTPUT
*/

#include "BoundaryType.H"
#include "ConvectionScheme.H"
#include "Vector.H"
#include "Mesh.H"
#include "Discretization.H"
#include "ScalarField.H"
#include "OutputData.H"
#include "OutputVTK.H"

#include <chrono>
#include <cstddef>
#include <string>
#include <iostream>
#include <Eigen/Dense>
#include <fstream>
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
        READ CONFIG FILE DATA INPUTS
    */
    Config cfg;
    cfg.readFile(argv[1]);
    const Setting& root = cfg.getRoot();

    // Set simulation data
    const Setting& simulationSetup = root["simulationSetup"];
    std::vector<int> nCells(3);                                                                       // No. of cells
    nCells = {simulationSetup.lookup("nCellsx"), simulationSetup.lookup("nCellsy"), simulationSetup.lookup("nCellsz")};
    int maxIterations = simulationSetup.lookup("maxIterations");                                   // No. of iterations
    ConvectionScheme fluxScheme = stringToConvectionScheme(simulationSetup.lookup("fluxScheme"));                                 // Flux scheme for convection discretization

    // Set initial data
    const Setting& inputParams = root["inputParams"];
    int Re = inputParams.lookup("Re");                                                          // Reynolds number
    double ux = inputParams.lookup("velX"), uy = inputParams.lookup("velY");                 // Initial velocity
    double alpha_P = inputParams.lookup("alpha_P"), alpha_u = inputParams.lookup("alpha_u"); // Under-relaxation factors
    double L = inputParams.lookup("L");                                                           // Cavity dimensions
    double nu = sqrt(ux * ux + uy * uy) * L / Re;                                             // Calculate nu, where, Re = u * L / nu

    // Read in file names for mesh generation
    const Setting& fileLocations = root["fileLocations"];
    Mesh cavity_mesh(fileLocations.lookup("points"), fileLocations.lookup("faces"), fileLocations.lookup("cells"), fileLocations.lookup("boundaries"));
    std::vector<Cell> myListOfAllCls = cavity_mesh.getListOfCells(); // Get list of cells in mesh
    std::vector<Face> myListOfAllFcs = cavity_mesh.getListOfFaces(); // Get list of all faces in the mesh
    std::vector<Point> myListOfAllPts = cavity_mesh.getListOfPoints(); // Get list of all points in the mesh

    // Get output location (to store data files)
    std::string outputFolder = fileLocations.lookup("outputFolder");
    

    /*
        SET UP FIELD VARIABLES AND THEIR BOUNDARY CONDITIONS 
    */
    // 1. Velocity field
    VectorField velocity(myListOfAllCls);
    velocity.setInitialCondition(simulationSetup.lookup("velXInit"), simulationSetup.lookup("velYInit"), simulationSetup.lookup("velZInit"));
    const Setting& velocityBCs = root["velocityBCs"];
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

    // 2. Pressure field
    ScalarField pressure(myListOfAllCls);  
    pressure.setInitialCondition(simulationSetup.lookup("pressureInit"));
    const Setting& pressureBCs = root["pressureBCs"];
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



    // /*
    //     FILE OUT FOR RESIDUAL NORM (need to uncomment another line in the code!)
    // */
    // // std::ofstream outputResNorm("gradedCavity/Results/ResNorm.txt");
    // // outputResNorm << "Iteration P_norm Vx_norm Vy_norm Vz_norm" << std::endl;



    // // /*
    // //     FILE OUT FOR CONVERGENCE DATA
    // // */
    // // // Set up output for velocity and pressure convergence data
    // // std::ofstream outputConvergenceData("40x40Cavity/Results/ConvergenceData.txt");
    // // outputConvergenceData << "Iteration VelXAfterMomPred VelYAfterMomPred PressureAfterPressurePred VelXAfterCorr VelYAfterCorr PressureAfterCorr" << std::endl;
    // // Eigen::VectorXd velXAfterMomPred, velYAfterMomPred, velXAfterCorr, velYAfterCorr, presAfterPresPred, presAfterCorr;



    /* 
            SIMPLE ALGORITHM LOOP
    */
    for (int iteration = 0; iteration < maxIterations; iteration++)
    {
        /*
            SET UP AND SOLVE THE MOMENTUM EQUATION
        */
        evaluateConvectionTerm(myListOfAllFcs, velocity, velocity, fluxScheme);
        evaluateDiffusionTerm(myListOfAllFcs, velocity, nu);

        // Under relax velocity
        underRelaxImplicit(velocity, alpha_u);

        // Solve the momentum equation 
        velocity.solveEquation();


        // // Output velocity data after momentum predictor - for convergence data
        // if (iteration > 0)
        // {
        //     double diffX = 0.0, diffY = 0.0;
        //     double max_val_X = 0.0, max_val_Y = 0.0;

        //     for (int i = 0; i < velocity.getMatrixSize(); i++)
        //     {
        //         diffX = fabs(velocity.variableX(i) - velXAfterMomPred(i));
        //         diffY = fabs(velocity.variableY(i) - velYAfterMomPred(i));

        //         if (i == 0)
        //         {
        //             max_val_X = diffX;
        //             max_val_Y = diffY;
        //         }

        //         if (diffX > max_val_X)
        //         {
        //             max_val_X = diffX;
        //         }

        //         if (diffY > diffY)
        //         {
        //             max_val_Y = diffY;
        //         }
        //     }
        //     outputConvergenceData << iteration << " " << max_val_X << " " << max_val_Y;
        // }
        // velXAfterMomPred = velocity.variableX;
        // velYAfterMomPred = velocity.variableY;

        // Reset diagonal of momentum coeff matrix
        for (int i = 0; i < velocity.getMatrixSize(); i++)
        {
            // velocity.coeffs(i, i) *= alpha_u;
            velocity.coeffs(i, i) = velocity.coeffs(i, i) * alpha_u;
        }



        /*
            SET UP AND SOLVE THE PRESSURE EQUATION
        */
        // Store p_old
        Eigen::VectorXd p_old = pressure.variable;  

        // Calculate F_pre and assign fluxes to faces
        for (size_t i = 0; i < myListOfAllFcs.size(); i++)
        {
            Face* f = &myListOfAllFcs[i];
            int idxP = f->getFaceOwner();
            int idxN = f->getFaceNeighbour();

            double F = getFlux(f, velocity, idxP, idxN);
            f->setFaceFlux(F);
        }

        // Set div u contribution to source matrix of p
        for (size_t i = 0; i < myListOfAllCls.size(); i++)
        {
            Cell c = myListOfAllCls[i];
            int cellIdx = c.getCellIndex();
            std::vector<Face *> myListOfFcs = c.getListOfFaces(); // Get faces for this cell
            double fluxForCell = 0.0;

            for (size_t j = 0; j < myListOfFcs.size(); j++)
            {
                Face* f = myListOfFcs[j]; // Pick a face from the list

                for (size_t m = 0; m < myListOfAllFcs.size(); m++)
                {
                    Face listFace = myListOfAllFcs[m];

                    if (f->getFaceIndex() == listFace.getFaceIndex())
                    {
                        if (listFace.getFaceOwner() == c.getCellIndex()) // face belongs to the cell
                        {
                            fluxForCell += listFace.getFaceFlux();
                        }
                        else // if the face belongs to a neighbour cell
                        {
                            fluxForCell += -listFace.getFaceFlux();
                        }
                    }
                }
            }

            pressure.source(cellIdx) += fluxForCell;
        }

        // Set up gamma for pressure field
        std::vector<double> gamma_p(myListOfAllCls.size());

        for (size_t i = 0; i < myListOfAllCls.size(); i++)
        {
            gamma_p[i] = 1.0 / velocity.coeffs(i, i);
        }

        // Evaluate diffusion term for pressure field and solve matrices
        evaluatePositiveDiffusionTerm(myListOfAllFcs, pressure, gamma_p); 

        // solve pressure eqn
        pressure.solveEquation();  


        // // Output pressure data after pressure predictor - for convergence data
        // if (iteration > 0)
        // {
        //     double diff = 0.0, max_val = 0.0;

        //     for (int i = 0; i < pressure.getMatrixSize(); i++)
        //     {
        //         diff = fabs(pressure.variable(i) - presAfterPresPred(i));

        //         if (i == 0)
        //         {
        //             max_val = diff;
        //         }

        //         if (diff > max_val)
        //         {
        //             max_val = diff;
        //         }
        //     }
        //     outputConvergenceData << " " << max_val;
        // }
        // presAfterPresPred = pressure.variable;



        /*
            CORRECT FACE FLUXES
        */
        // Vector to store corrected fluxes
        std::vector<double> F_corr(myListOfAllFcs.size());

        // Calculate F_corr and assign correct face fluxes
        for (size_t i = 0; i < myListOfAllFcs.size(); i++)
        {
            Face f = myListOfAllFcs[i];
            int idxP = f.getFaceOwner();
            int idxN = f.getFaceNeighbour();
            double F_pre = f.getFaceFlux();
            double a_P_interpolated = 0.0;
            double fx = f.getFaceInterpolationFactor();

            // Calc corrected flux -- assumes orthogonal mesh
            if (f.getIfBoundaryFace() == true)  // If face is a boundary
            {
                F_corr[i] = F_pre; // Since boundary condition on all boundaries is grad(P) = 0 
            }
            else  // if face is not a boundary
            {
                a_P_interpolated = (fx * (1.0 / velocity.coeffs(idxP, idxP))) + ((1.0 - fx) * (1.0 / velocity.coeffs(idxN, idxN)));
                F_corr[i] = F_pre - (a_P_interpolated * f.getFaceArea() * (pressure.variable(idxN) - pressure.variable(idxP)) / f.getFaceDelta());
            }
            
            // Update old face flux
            f.setFaceFlux(F_corr[i]);
        }



        /*
            CORRECT CELL CENTRE VELOCITIES
        */
        for (size_t i = 0; i < myListOfAllCls.size(); i++)
        {
            Cell c = myListOfAllCls[i];
            int idx = c.getCellIndex();
            Vector grad_P = getGradientGauss(pressure, idx, myListOfAllCls);

            velocity.variableX(idx) = velocity.variableX(idx) - ((1.0 / velocity.coeffs(idx, idx)) * grad_P[0] / c.getCellVolume());
            velocity.variableY(idx) = velocity.variableY(idx) - ((1.0 / velocity.coeffs(idx, idx)) * grad_P[1] / c.getCellVolume());
            velocity.variableZ(idx) = velocity.variableZ(idx) - ((1.0 / velocity.coeffs(idx, idx)) * grad_P[2] / c.getCellVolume());
        }

        // // Output velocity data after correction - for convergence data
        // if (iteration > 0)
        // {
        //     double diffX = 0.0, diffY = 0.0;
        //     double max_val_X = 0.0, max_val_Y = 0.0;

        //     for (int i = 0; i < velocity.getMatrixSize(); i++)
        //     {
        //         diffX = fabs(velocity.variableX(i) - velXAfterCorr(i));
        //         diffY = fabs(velocity.variableY(i) - velYAfterCorr(i));

        //         if (i == 0)
        //         {
        //             max_val_X = diffX;
        //             max_val_Y = diffY;
        //         }

        //         if (diffX > max_val_X)
        //         {
        //             max_val_X = diffX;
        //         }

        //         if (diffY > diffY)
        //         {
        //             max_val_Y = diffY;
        //         }
        //     }
        //     outputConvergenceData << " " << max_val_X << " " << max_val_Y;
        // }
        // velXAfterCorr = velocity.variableX;
        // velYAfterCorr = velocity.variableY;



        /*
            UNDER-RELAX PRESSURE
        */
        underRelaxExplicit(pressure, p_old, alpha_P);

        // // Output pressure data after pressure correction - for convergence data
        // if (iteration > 0)
        // {
        //     double diff = 0.0, max_val = 0.0;

        //     for (int i = 0; i < pressure.getMatrixSize(); i++)
        //     {
        //         diff = fabs(pressure.variable(i) - presAfterCorr(i));

        //         if (i == 0)
        //         {
        //             max_val = diff;
        //         }

        //         if (diff > max_val)
        //         {
        //             max_val = diff;
        //         }
        //     }

        //     outputConvergenceData << " " << max_val << std::endl;
        // }
        // presAfterCorr = pressure.variable;


        
        /*
            OUTPUT DATA REQUIRED FOR PLOTS FOR THIS ITERATION
        */
        // outputVariableDataForPlot(nCellsx, nCellsy, "40x40Cavity/Results/velocity", velocity, iteration);
        // outputVariableDataForPlot(nCellsx, nCellsy, "40x40Cavity/Results/pressure", pressure, iteration);
        // outputVectorMagnitudeForPlot(nCellsx, nCellsy, "40x40Cavity/Results/velocity", velocity, iteration);
        // outputNormalizedVectorCompForPlot(nCellsx, nCellsy, "40x40Cavity/Results/velocity", velocity, iteration);


        // Calculate and output residual data
        pressure.calculateResidual();
        pressure.calculateResidualNorm();
        velocity.calculateResidual();
        velocity.calculateResidualNorm();

        // outputAllResiduals(pressure, velocity, "40x40Cavity/Results/res", iteration);
        
        // // Calculate residual norms and output to file - needed for file output (see lines 117-118)
        // std::vector<double> velocityResNorm = velocity.getResidualNorm();
        // double pressureResNorm = pressure.getResidualNorm();
        // outputResNorm << iteration << " " << pressureResNorm << " " << velocityResNorm[0] << " " << velocityResNorm[1] << " " << velocityResNorm[2] << std::endl;


        // // Create VTK Output file
        // outputVTK_CavityData(myListOfAllPts, pressure, velocity, nCells, iteration, "40x40Cavity/VTK/");
        outputVTK_CavityData(myListOfAllPts, pressure, velocity, nCells, iteration, outputFolder);

        /*
            RESET COEFF AND SOURCE MATRICES BEFORE NEXT TIME STEP
        */
        velocity.resetMatrices();
        pressure.resetMatrices();

        // Print iteration number at the end of the loop
        std::cout << "Iteration " << iteration << std::endl;
    }

    // // Output res norm data
    // outputAllResidualNorms(pressure, velocity, "40x40Cavity/Results/");

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