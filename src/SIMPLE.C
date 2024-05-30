/*
    SIMPLE ALGORITHM
*/

#include "SIMPLE.H"

#include "ConvectionScheme.H"
#include "Vector.H"
#include "MRFMesh.H"
#include "SRFMesh.H"
#include "MovingMeshDiscretization.H"
#include "VectorField.H"
#include "ScalarField.H"
#include "OutputVTK.H"
#include "OutputData.H"

#include <iostream>
#include <cstddef>
#include <Eigen/Dense>
#include <libconfig.h++>
#include <fstream>

using namespace libconfig;


void evaluateSIMPLE(MRFMesh& mesh, VectorField& velocity, ScalarField& pressure, const Setting& root)
{   
    // Get mesh related info
    std::vector<Cell> myListOfAllCls = mesh.getListOfCells(); // Get list of cells in mesh
    std::vector<Face> myListOfAllFcs = mesh.getListOfFaces(); // Get list of all faces in the mesh
    std::vector<Point> myListOfAllPts = mesh.getListOfPoints(); // Get list of all points in the mesh

    // Get config file data
    const Setting& inputParams = root["inputParams"];
    const Setting& simulationSetup = root["simulationSetup"];
    const Setting& fileLocations = root["fileLocations"];

    // Get simulationData info
    const Setting& simulationData = root["simulationData"];
    std::string outputDataFolder = simulationData.lookup("outputFolder"); 

    // By default, no output files are generated (apart from VTK)
    bool outputAllData = false, outputConvData = false, outputResiduals = false, outputResidualNorm = false;
    if (simulationData.exists("outputAllData"))
    {
        outputAllData = simulationData.lookup("outputAllData");
    }
    else
    {
        outputConvData = simulationData.lookup("outputConvergenceData");
        outputResiduals = simulationData.lookup("outputResiduals");
        outputResidualNorm = simulationData.lookup("outputResidualNorm");
    }

    bool isMRF = simulationSetup.lookup("isMRF");
    std::string outputFolder = fileLocations.lookup("outputFolder");                               // Get output location (to store data files)
    int maxIterations = simulationSetup.lookup("maxIterations");                                   // No. of iterations
    ConvectionScheme fluxScheme = stringToConvectionScheme(simulationSetup.lookup("fluxScheme"));                                 // Flux scheme for convection discretization               
    double alpha_P = inputParams.lookup("alpha_P"), alpha_u = inputParams.lookup("alpha_u"); // Under-relaxation factors  
    std::string movingZone = simulationSetup.lookup("movingZone");                                           

    // Get or Calculate nu
    double nu = 0.0;
    if (inputParams.exists("nu"))
    {
        nu = inputParams.lookup("nu");
    }
    else if (inputParams.exists("Re") && inputParams.exists("inputVelocity") && inputParams.exists("L"))
    {
        int Re = inputParams.lookup("Re");  // Reynolds number
        double L = inputParams.lookup("L"); // Cavity dimensions

        std::vector<double> u;                    // Input velocity
        const Setting& inputVelSett = simulationSetup["inputVelocity"];
        for (int i = 0; i < 3; i++)
        {
            u.push_back(inputVelSett[i]);
        }

        // Calculate nu, where, Re = u * L / nu
        nu = sqrt(u[0] * u[0] + u[1] * u[1] + u[2] * u[2]) * L / Re; 
    }
    else
    {
        std::cerr << "Insufficient data provided! Cannot calculate nu." << std::endl;
        exit(1);
    }

    // Get omega and origin from config file
    Vector omega, origin;
    std::vector<double> omegaVect, originVect;
    const Setting& omegaSett = simulationSetup["omega"];
    const Setting& orginSett = simulationSetup["origin"];

    for (int i = 0; i < 3; i++)
    {
        omegaVect.push_back(omegaSett[i]);
        originVect.push_back(orginSett[i]);
    }
    omega = {omegaVect[0], omegaVect[1], omegaVect[2]};
    origin = {originVect[0], originVect[1], originVect[2]};


    // Create VTK Output file for initial condition (i.e., iteration 0)
    outputVTK_MRFData(mesh, pressure, velocity, 0, outputFolder);



    /*
        FILE OUT FOR CONVERGENCE DATA - Set up output for velocity and pressure convergence data
    */
    std::ofstream outputConvergenceData(outputDataFolder + "ConvergenceData.txt");
    if (outputAllData == true || outputConvData == true)
    {
        outputConvergenceData << "Iteration VelXAfterMomPred VelYAfterMomPred PressureAfterPressurePred VelXAfterCorr VelYAfterCorr PressureAfterCorr" << std::endl;
    }
    Eigen::VectorXd velXAfterMomPred, velYAfterMomPred, velXAfterCorr, velYAfterCorr, presAfterPresPred, presAfterCorr;


    /*
        SIMPLE ALGORITHM
    */
    for (int iteration = 1; iteration <= maxIterations; iteration++)
    {
        /*
            SET UP AND SOLVE THE MOMENTUM EQUATION
        */
        evaluateConvectionTerm(myListOfAllFcs, velocity, velocity, fluxScheme, omega, origin);
        evaluateDiffusionTerm(myListOfAllFcs, velocity, nu, omega, origin);

        // Add source term contribution for MRF
        for (int i = 0; i < velocity.getMatrixSize(); i++)
        {
            Vector vel = {velocity.variableX(i), velocity.variableY(i), velocity.variableZ(i)};
            Vector MRF_source = omega.crossMult(vel) * -1.0 * myListOfAllCls[i].getCellVolume();

            velocity.addSourceContribution(i, MRF_source[0], MRF_source[1], MRF_source[2]);
        }

        // Under relax velocity
        underRelaxImplicit(velocity, alpha_u);

        // Solve the momentum equation 
        velocity.solveEquation();


        // Output velocity data after momentum predictor - for convergence data
        if (outputAllData == true || outputConvData == true)
        {
            if (iteration > 1)
            {
                double diffX = 0.0, diffY = 0.0;
                double max_val_X = 0.0, max_val_Y = 0.0;

                for (int i = 0; i < velocity.getMatrixSize(); i++)
                {
                    diffX = fabs(velocity.variableX(i) - velXAfterMomPred(i));
                    diffY = fabs(velocity.variableY(i) - velYAfterMomPred(i));

                    if (i == 0)
                    {
                        max_val_X = diffX;
                        max_val_Y = diffY;
                    }

                    if (diffX > max_val_X)
                    {
                        max_val_X = diffX;
                    }

                    if (diffY > diffY)
                    {
                        max_val_Y = diffY;
                    }
                }
                outputConvergenceData << iteration << " " << max_val_X << " " << max_val_Y;
            }
            velXAfterMomPred = velocity.variableX;
            velYAfterMomPred = velocity.variableY;
        }


        // Reset diagonal of momentum coeff matrix
        for (int i = 0; i < velocity.getMatrixSize(); i++)
        {
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

            double F = getFlux(f, myListOfAllCls, velocity, idxP, idxN, omega, origin, isMRF);
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


        // Output pressure data after pressure predictor - for convergence data
        if (outputAllData == true || outputConvData == true)
        {
            if (iteration > 1)
            {
                double diff = 0.0, max_val = 0.0;

                for (int i = 0; i < pressure.getMatrixSize(); i++)
                {
                    diff = fabs(pressure.variable(i) - presAfterPresPred(i));

                    if (i == 0)
                    {
                        max_val = diff;
                    }

                    if (diff > max_val)
                    {
                        max_val = diff;
                    }
                }
                outputConvergenceData << " " << max_val;
            }
            presAfterPresPred = pressure.variable;
        }


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


        // Output velocity data after correction - for convergence data
        if (outputAllData == true || outputConvData == true)
        {
            if (iteration > 1)
            {
                double diffX = 0.0, diffY = 0.0;
                double max_val_X = 0.0, max_val_Y = 0.0;

                for (int i = 0; i < velocity.getMatrixSize(); i++)
                {
                    diffX = fabs(velocity.variableX(i) - velXAfterCorr(i));
                    diffY = fabs(velocity.variableY(i) - velYAfterCorr(i));

                    if (i == 0)
                    {
                        max_val_X = diffX;
                        max_val_Y = diffY;
                    }

                    if (diffX > max_val_X)
                    {
                        max_val_X = diffX;
                    }

                    if (diffY > diffY)
                    {
                        max_val_Y = diffY;
                    }
                }
                outputConvergenceData << " " << max_val_X << " " << max_val_Y;
            }
            velXAfterCorr = velocity.variableX;
            velYAfterCorr = velocity.variableY;
        }


        /*
            UNDER-RELAX PRESSURE
        */
        underRelaxExplicit(pressure, p_old, alpha_P);

        // Output pressure data after pressure correction - for convergence data
        if (outputAllData == true || outputConvData == true)
        {
            if (iteration > 1)
            {
                double diff = 0.0, max_val = 0.0;

                for (int i = 0; i < pressure.getMatrixSize(); i++)
                {
                    diff = fabs(pressure.variable(i) - presAfterCorr(i));

                    if (i == 0)
                    {
                        max_val = diff;
                    }

                    if (diff > max_val)
                    {
                        max_val = diff;
                    }
                }

                outputConvergenceData << " " << max_val << std::endl;
            }
            presAfterCorr = pressure.variable;
        }

        // Calculate and output residual data
        pressure.calculateResidual();
        pressure.calculateResidualNorm();
        velocity.calculateResidual();
        velocity.calculateResidualNorm();

        // Output residuals for all fields
        if (outputAllData == true || outputResiduals == true)
        {
            outputAllResiduals(pressure, velocity, outputDataFolder + "res", iteration);
        }

        // Create VTK Output file
        outputVTK_MRFData(mesh, pressure, velocity, iteration, outputFolder);


        /*
            RESET COEFF AND SOURCE MATRICES BEFORE NEXT TIME STEP
        */
        velocity.resetMatrices();
        pressure.resetMatrices();

        // Print iteration number at the end of the loop
        std::cout << "Iteration " << iteration << std::endl;
    }

    // Output res norm data
    if (outputAllData == true || outputResidualNorm == true)
    {
        outputAllResidualNorms(pressure, velocity, outputDataFolder);
    }
}



void evaluateSIMPLE(SRFMesh& mesh, VectorField& velocity, ScalarField& pressure, const Setting& root)
{   
    // Get mesh related info
    std::vector<Cell> myListOfAllCls = mesh.getListOfCells(); // Get list of cells in mesh
    std::vector<Face> myListOfAllFcs = mesh.getListOfFaces(); // Get list of all faces in the mesh
    std::vector<Point> myListOfAllPts = mesh.getListOfPoints(); // Get list of all points in the mesh

    // Get config file data
    const Setting& inputParams = root["inputParams"];
    const Setting& simulationSetup = root["simulationSetup"];
    const Setting& fileLocations = root["fileLocations"];

    // Get simulationData info
    const Setting& simulationData = root["simulationData"];
    std::string outputDataFolder = simulationData.lookup("outputFolder"); 

    // By default, no output files are generated (apart from VTK)
    bool outputAllData = false, outputConvData = false, outputResiduals = false, outputResidualNorm = false;
    if (simulationData.exists("outputAllData"))
    {
        outputAllData = simulationData.lookup("outputAllData");
    }
    else
    {
        outputConvData = simulationData.lookup("outputConvergenceData");
        outputResiduals = simulationData.lookup("outputResiduals");
        outputResidualNorm = simulationData.lookup("outputResidualNorm");
    }

    bool isMRF = simulationSetup.lookup("isMRF");
    std::string outputFolder = fileLocations.lookup("outputFolder");                               // Get output location (to store data files)
    int maxIterations = simulationSetup.lookup("maxIterations");                                   // No. of iterations
    ConvectionScheme fluxScheme = stringToConvectionScheme(simulationSetup.lookup("fluxScheme"));                                 // Flux scheme for convection discretization               
    double alpha_P = inputParams.lookup("alpha_P"), alpha_u = inputParams.lookup("alpha_u"); // Under-relaxation factors                                             

    // Get or Calculate nu
    double nu = 0.0;
    if (inputParams.exists("nu"))
    {
        nu = inputParams.lookup("nu");
    }
    else if (inputParams.exists("Re") && inputParams.exists("inputVelocity") && inputParams.exists("L"))
    {
        int Re = inputParams.lookup("Re");  // Reynolds number
        double L = inputParams.lookup("L"); // Cavity dimensions

        std::vector<double> u;                    // Input velocity
        const Setting& inputVelSett = simulationSetup["inputVelocity"];
        for (int i = 0; i < 3; i++)
        {
            u.push_back(inputVelSett[i]);
        }

        // Calculate nu, where, Re = u * L / nu
        nu = sqrt(u[0] * u[0] + u[1] * u[1] + u[2] * u[2]) * L / Re; 
    }
    else
    {
        std::cerr << "Insufficient data provided! Cannot calculate nu." << std::endl;
        exit(1);
    }

    // Get omega and origin from config file
    Vector omega, origin;
    std::vector<double> omegaVect, originVect;
    const Setting& omegaSett = simulationSetup["omega"];
    const Setting& orginSett = simulationSetup["origin"];

    for (int i = 0; i < 3; i++)
    {
        omegaVect.push_back(omegaSett[i]);
        originVect.push_back(orginSett[i]);
    }
    omega = {omegaVect[0], omegaVect[1], omegaVect[2]};
    origin = {originVect[0], originVect[1], originVect[2]};


    // Create VTK Output file for initial condition (i.e., iteration 0)
    outputVTK_SRFData(mesh, pressure, velocity, 0, outputFolder);



    /*
        FILE OUT FOR CONVERGENCE DATA - Set up output for velocity and pressure convergence data
    */
    std::ofstream outputConvergenceData(outputDataFolder + "ConvergenceData.txt");
    if (outputAllData == true || outputConvData == true)
    {
        outputConvergenceData << "Iteration VelXAfterMomPred VelYAfterMomPred PressureAfterPressurePred VelXAfterCorr VelYAfterCorr PressureAfterCorr" << std::endl;
    }
    Eigen::VectorXd velXAfterMomPred, velYAfterMomPred, velXAfterCorr, velYAfterCorr, presAfterPresPred, presAfterCorr;


    /*
        SIMPLE ALGORITHM
    */
    for (int iteration = 1; iteration <= maxIterations; iteration++)
    {
        /*
            SET UP AND SOLVE THE MOMENTUM EQUATION
        */
        // evaluateTimeDerivative(myListOfAllCls, velocity, 0.1);
        evaluateConvectionTerm(myListOfAllFcs, velocity, velocity, fluxScheme, omega, origin);
        evaluateDiffusionTerm(myListOfAllFcs, velocity, nu, omega, origin);

        // Add source term contributions - for coriolis force
        for (int i = 0; i < velocity.getMatrixSize(); i++)
        {
            Vector vel = {velocity.variableX(i), velocity.variableY(i), velocity.variableZ(i)};
            Vector coriolis = omega.crossMult(vel) * -2.0 * myListOfAllCls[i].getCellVolume();

            velocity.addSourceContribution(i, coriolis[0], coriolis[1], coriolis[2]);
        }

        // Add source term contributions - for centrifugal force
        for (int i = 0; i < velocity.getMatrixSize(); i++)
        {
            Vector centrifugal = omega.crossMult(omega.crossMult(myListOfAllCls[i].getCellRadius())) * -1.0 * myListOfAllCls[i].getCellVolume();

            velocity.addSourceContribution(i, centrifugal[0], centrifugal[1], centrifugal[2]);
        }

        // Under relax velocity
        underRelaxImplicit(velocity, alpha_u);

        // Solve the momentum equation 
        velocity.solveEquation();


        // Output velocity data after momentum predictor - for convergence data
        if (outputAllData == true || outputConvData == true)
        {
            if (iteration > 1)
            {
                double diffX = 0.0, diffY = 0.0;
                double max_val_X = 0.0, max_val_Y = 0.0;

                for (int i = 0; i < velocity.getMatrixSize(); i++)
                {
                    diffX = fabs(velocity.variableX(i) - velXAfterMomPred(i));
                    diffY = fabs(velocity.variableY(i) - velYAfterMomPred(i));

                    if (i == 0)
                    {
                        max_val_X = diffX;
                        max_val_Y = diffY;
                    }

                    if (diffX > max_val_X)
                    {
                        max_val_X = diffX;
                    }

                    if (diffY > diffY)
                    {
                        max_val_Y = diffY;
                    }
                }
                outputConvergenceData << iteration << " " << max_val_X << " " << max_val_Y;
            }
            velXAfterMomPred = velocity.variableX;
            velYAfterMomPred = velocity.variableY;
        }


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

            double F = getFlux(f, myListOfAllCls, velocity, idxP, idxN, omega, origin, isMRF);
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


        // Output pressure data after pressure predictor - for convergence data
        if (outputAllData == true || outputConvData == true)
        {
            if (iteration > 1)
            {
                double diff = 0.0, max_val = 0.0;

                for (int i = 0; i < pressure.getMatrixSize(); i++)
                {
                    diff = fabs(pressure.variable(i) - presAfterPresPred(i));

                    if (i == 0)
                    {
                        max_val = diff;
                    }

                    if (diff > max_val)
                    {
                        max_val = diff;
                    }
                }
                outputConvergenceData << " " << max_val;
            }
            presAfterPresPred = pressure.variable;
        }


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


        // Output velocity data after correction - for convergence data
        if (outputAllData == true || outputConvData == true)
        {
            if (iteration > 1)
            {
                double diffX = 0.0, diffY = 0.0;
                double max_val_X = 0.0, max_val_Y = 0.0;

                for (int i = 0; i < velocity.getMatrixSize(); i++)
                {
                    diffX = fabs(velocity.variableX(i) - velXAfterCorr(i));
                    diffY = fabs(velocity.variableY(i) - velYAfterCorr(i));

                    if (i == 0)
                    {
                        max_val_X = diffX;
                        max_val_Y = diffY;
                    }

                    if (diffX > max_val_X)
                    {
                        max_val_X = diffX;
                    }

                    if (diffY > diffY)
                    {
                        max_val_Y = diffY;
                    }
                }
                outputConvergenceData << " " << max_val_X << " " << max_val_Y;
            }
            velXAfterCorr = velocity.variableX;
            velYAfterCorr = velocity.variableY;
        }


        /*
            UNDER-RELAX PRESSURE
        */
        underRelaxExplicit(pressure, p_old, alpha_P);

        // Output pressure data after pressure correction - for convergence data
        if (outputAllData == true || outputConvData == true)
        {
            if (iteration > 1)
            {
                double diff = 0.0, max_val = 0.0;

                for (int i = 0; i < pressure.getMatrixSize(); i++)
                {
                    diff = fabs(pressure.variable(i) - presAfterCorr(i));

                    if (i == 0)
                    {
                        max_val = diff;
                    }

                    if (diff > max_val)
                    {
                        max_val = diff;
                    }
                }

                outputConvergenceData << " " << max_val << std::endl;
            }
            presAfterCorr = pressure.variable;
        }

        // Calculate and output residual data
        pressure.calculateResidual();
        pressure.calculateResidualNorm();
        velocity.calculateResidual();
        velocity.calculateResidualNorm();

        // Output residuals for all fields
        if (outputAllData == true || outputResiduals == true)
        {
            outputAllResiduals(pressure, velocity, outputDataFolder + "res", iteration);
        }

        // Create VTK Output file
        outputVTK_SRFData(mesh, pressure, velocity, iteration, outputFolder);


        /*
            RESET COEFF AND SOURCE MATRICES BEFORE NEXT TIME STEP
        */
        velocity.resetMatrices();
        pressure.resetMatrices();

        // Print iteration number at the end of the loop
        std::cout << "Iteration " << iteration << std::endl;
    }

    // Output res norm data
    if (outputAllData == true || outputResidualNorm == true)
    {
        outputAllResidualNorms(pressure, velocity, outputDataFolder);
    }
}