#include "OutputVTK.H"

#include "ScalarField.H"
#include "VectorField.H"
#include "Vector.H"
#include "Point.H"
#include "Face.H"
#include "Cell.H"
#include "Mesh.H"
#include "MRFMesh.H"
#include "SRFMesh.H"

#include <cstddef>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <array>

#include <vtkCellType.h>
#include <vtkGenericCell.h>
#include <vtkIdList.h>
#include <vtkSmartPointer.h>
#include <vtkPoints.h>
#include <vtkCellData.h>
#include <vtkCellArray.h>
#include <vtkType.h>
#include <vtkXMLStructuredGridWriter.h>
#include <vtkStructuredGrid.h>
#include <vtkHexahedron.h>
#include <vtkDoubleArray.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPolygon.h>
#include <vtkXMLWriter.h>
#include <vtkPolyData.h>
#include <vtkXMLPolyDataWriter.h>


// VTK output for lid-driven cavity
void outputVTK_CavityData(std::vector<Point>& listOfAllPoints, ScalarField& pressure, VectorField& velocity, std::vector<int> nCells, int iteration, std::string filename)
{
    int matrixSize = pressure.getMatrixSize();  // Is the same as total number of cells in the domain

    /*
        CREATE MESH
    */
    vtkSmartPointer<vtkStructuredGrid> grid = vtkSmartPointer<vtkStructuredGrid>::New();
    grid->SetDimensions(nCells[0]+1, nCells[1]+1, nCells[2]+1);

    // Create points 
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

    // Add each point from the list of points in the mesh
    for (size_t i = 0; i < listOfAllPoints.size(); i++)
    {
        Vector p = listOfAllPoints[i].pointToVector();
        points->InsertNextPoint(p[0], p[1], p[2]);
    }

    // Add points to the grid
    grid->SetPoints(points);


    /*
        ADD SCALAR DATA
    */
    // Define names for all scalar variables
    std::vector<std::string> scalarFieldNames = {"Pressure", "Velocity_X", "Velocity_Y", "Velocity_Z"};

    // Vector of vectors to store all scalar data for all 4 variables -> pressure, velX, velY, velZ
    std::vector<std::vector<double>> scalarValues(4, std::vector<double> (matrixSize));

    // Store data in scalarValues
    for (int i = 0; i < matrixSize; i++)
    {
        scalarValues[0][i] = pressure.variable(i);
        scalarValues[1][i] = velocity.variableX(i);
        scalarValues[2][i] = velocity.variableY(i);
        scalarValues[3][i] = velocity.variableZ(i);
    }

    // For every scalar variable, set cell centre values 
    for (int var = 0; var < 4; var++)
    {
        // Create a VTK data array for the scalar variable
        vtkSmartPointer<vtkDoubleArray> scalarData = vtkSmartPointer<vtkDoubleArray>::New();

        // Set the name of the variable
        scalarData->SetName(scalarFieldNames[var].c_str()); 

        // Set values for each cell
        for (int i = 0; i < matrixSize; i++)
        {
            scalarData->InsertNextValue(scalarValues[var][i]);
        }

        grid->GetCellData()->AddArray(scalarData);
    }


    /*
        ADD VECTOR DATA
    */
    vtkSmartPointer<vtkDoubleArray> velocityData = vtkSmartPointer<vtkDoubleArray>::New();
    velocityData->SetName("Velocity");
    velocityData->SetNumberOfComponents(3); // Assumes 3 vector components

    for (int i = 0; i < velocity.getMatrixSize(); i++)
    {
        std::array<double, 3> velData = {velocity.variableX(i), velocity.variableY(i), velocity.variableZ(i)};
        velocityData->InsertNextTuple(velData.data());
    }

    grid->GetCellData()->SetVectors(velocityData);


    /*
        OUTPUT DATA TO VTK XML FILE
    */
    vtkSmartPointer<vtkXMLStructuredGridWriter> writer = vtkSmartPointer<vtkXMLStructuredGridWriter>::New();
    std::stringstream ss;
    ss << filename + "Cavity_Itr_" << iteration << ".vts";
    std::string fileName = ss.str();
    writer->SetFileName(fileName.c_str()); 
    writer->SetDataModeToBinary();
    writer->SetInputData(grid);      
    writer->Write(); 
}



// VTK output for MRF data 
// void outputVTK_MRFData(std::vector<Point>& listOfAllPoints, std::vector<Face>& listOfAllFaces, std::vector<Cell>& listOfAllCells, ScalarField& pressure, VectorField& velocity, int iteration, std::string filename)
void outputVTK_MRFData(MRFMesh& mesh, ScalarField& pressure, VectorField& velocity, int iteration, std::string filename)
{
    // Get required mesh details
    std::vector<Point> listOfAllPoints = mesh.getListOfPoints();
    std::vector<Cell> listOfAllCells = mesh.getListOfCells();

    int matrixSize = pressure.getMatrixSize();  // Is the same as total number of cells in the domain


    // Create an unstructured VTK grid
    vtkSmartPointer<vtkUnstructuredGrid> grid = vtkSmartPointer<vtkUnstructuredGrid>::New();

    // Create points 
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

    // Add each point from the list of points in the mesh
    for (size_t i = 0; i < listOfAllPoints.size(); i++)
    {
        Vector p = listOfAllPoints[i].pointToVector();
        points->InsertNextPoint(p[0], p[1], p[2]);
    }

    // Add points to the grid
    grid->SetPoints(points);



    // Create a list for all cells in the grid
    vtkSmartPointer<vtkCellArray> allCells = vtkSmartPointer<vtkCellArray>::New(); // list of all cells

    // Add cells to the grid
    for (size_t i = 0; i < listOfAllCells.size(); i++)
    {
        // Get a cell
        Cell c = listOfAllCells[i];

        // Get list of faces for this cell
        std::vector<Face *> faces = c.getListOfFaces();

        // Vector to store points of a cell
        std::vector<Point> pts;

        // Vector to store indices of points
        std::vector<int> ptIdx;

        // Get points of cell
        for (size_t j = 0; j < faces.size(); j++)
        {
            // Get a face
            Face* f = faces[j];      
            pts = f->getListOfPoints();

            // Store points in vector of pts for this cell
            if (j == 0)
            {
                for (size_t k = 0; k < pts.size(); k++)
                {
                    ptIdx.push_back(pts[k].getPointIndex());
                }
            }
            else 
            {
                bool found;

                for (size_t k = 0; k < f->getListOfPoints().size(); k++)
                {
                    Point p = f->getListOfPoints()[k];
                    int idx = p.getPointIndex();

                    found = (std::find(ptIdx.begin(), ptIdx.end(), idx) != ptIdx.end());

                    if (found == true)
                    {
                        break;
                    }
                }

                if (found == false)
                {
                    for (size_t k = 0; k < pts.size(); k++)
                    {
                        ptIdx.push_back(pts[k].getPointIndex());
                    }
                }
            }
        }



        // Add cell to the grid - assumes tetrahedral cells
        vtkSmartPointer<vtkHexahedron> cell = vtkSmartPointer<vtkHexahedron>::New();

        for (size_t j = 0; j < 8; j++)
        {
            cell->GetPointIds()->SetId(j, ptIdx[j]);
        }

        allCells->InsertNextCell(cell);
    }

    grid->SetCells(VTK_HEXAHEDRON, allCells);



    /*
        ADD SCALAR DATA
    */
    // Define names for all scalar variables
    std::vector<std::string> scalarFieldNames = {"Pressure", "Velocity_X", "Velocity_Y", "Velocity_Z"};

    // Vector of vectors to store all scalar data for all 4 variables -> pressure, velX, velY, velZ
    std::vector<std::vector<double>> scalarValues(4, std::vector<double> (matrixSize));

    // Store data in scalarValues
    for (int i = 0; i < matrixSize; i++)
    {
        scalarValues[0][i] = pressure.variable(i) / listOfAllCells[i].getCellVolume(); // For pressure, divide by volume for output;
        scalarValues[1][i] = velocity.variableX(i);
        scalarValues[2][i] = velocity.variableY(i);
        scalarValues[3][i] = velocity.variableZ(i);
    }

    // For every scalar variable, set cell centre values 
    for (int var = 0; var < 4; var++)
    {
        // Create a VTK data array for the scalar variable
        vtkSmartPointer<vtkDoubleArray> scalarData = vtkSmartPointer<vtkDoubleArray>::New();

        // Set the name of the variable
        scalarData->SetName(scalarFieldNames[var].c_str()); 

        // Set values for each cell
        for (int i = 0; i < matrixSize; i++)
        {
            scalarData->InsertNextValue(scalarValues[var][i]);
        }

        grid->GetCellData()->AddArray(scalarData);
    }


    /*
        ADD VECTOR DATA
    */
    vtkSmartPointer<vtkDoubleArray> velocityData = vtkSmartPointer<vtkDoubleArray>::New();
    velocityData->SetName("Velocity");
    velocityData->SetNumberOfComponents(3); // Assumes 3 vector components

    for (int i = 0; i < velocity.getMatrixSize(); i++)
    {
        std::array<double, 3> velData = {velocity.variableX(i), velocity.variableY(i), velocity.variableZ(i)};
        velocityData->InsertNextTuple(velData.data());
    }

    grid->GetCellData()->SetVectors(velocityData);



    /*
        Output VTK XML file
    */
    // Write the grid to a VTK XML unstructured grid file
    vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();

    std::stringstream ss;
    ss << filename + "MRF_Itr_" << iteration << ".vtu";
    std::string fileName = ss.str();
    writer->SetFileName(fileName.c_str()); 

    writer->SetDataModeToBinary();
    writer->SetInputData(grid);
    writer->Write(); 
}



// VTK output for SRF data 
void outputVTK_SRFData(SRFMesh& mesh, ScalarField& pressure, VectorField& velocity, int iteration, std::string filename)
// void outputVTK_SRFData(SRFMesh& mesh, int iteration, std::string filename)
{
    // Get required mesh details
    std::vector<Point> listOfAllPoints = mesh.getListOfPoints();
    std::vector<Cell> listOfAllCells = mesh.getListOfCells();

    int matrixSize = pressure.getMatrixSize();  // Is the same as total number of cells in the domain

    // Create an unstructured VTK grid
    vtkSmartPointer<vtkUnstructuredGrid> grid = vtkSmartPointer<vtkUnstructuredGrid>::New();

    // Create points 
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

    // Add each point from the list of points in the mesh
    for (size_t i = 0; i < listOfAllPoints.size(); i++)
    {
        Vector p = listOfAllPoints[i].pointToVector();
        points->InsertNextPoint(p[0], p[1], p[2]);
    }

    // Add points to the grid
    grid->SetPoints(points);


    // Create a list for all cells in the grid
    vtkSmartPointer<vtkCellArray> allCells = vtkSmartPointer<vtkCellArray>::New(); // list of all cells

    // Add cells to the grid
    for (size_t i = 0; i < listOfAllCells.size(); i++)
    {
        // Get a cell
        Cell c = listOfAllCells[i];

        // Get list of faces for this cell
        std::vector<Face *> faces = c.getListOfFaces();


        // Vector to store points of a cell
        std::vector<Point> pts;

        // Vector to store indices of points
        std::vector<int> ptIdx;

        // Get points of cell
        for (size_t j = 0; j < faces.size(); j++)
        {
            // Get a face
            Face* f = faces[j];      
            pts = f->getListOfPoints();

            // Store points in vector of pts for this cell
            if (j == 0)
            {
                for (size_t k = 0; k < pts.size(); k++)
                {
                    ptIdx.push_back(pts[k].getPointIndex());
                }
            }
            else 
            {
                bool found;

                for (size_t k = 0; k < f->getListOfPoints().size(); k++)
                {
                    Point p = f->getListOfPoints()[k];
                    int idx = p.getPointIndex();

                    found = (std::find(ptIdx.begin(), ptIdx.end(), idx) != ptIdx.end());

                    if (found == true)
                    {
                        break;
                    }
                }

                if (found == false)
                {
                    for (size_t k = 0; k < pts.size(); k++)
                    {
                        ptIdx.push_back(pts[k].getPointIndex());
                    }
                }
            }
        }

        vtkSmartPointer<vtkHexahedron> cell = vtkSmartPointer<vtkHexahedron>::New();

        for (size_t j = 0; j < 8; j++)
        {
            cell->GetPointIds()->SetId(j, ptIdx[j]);
        }

        allCells->InsertNextCell(cell);
    }

    grid->SetCells(VTK_HEXAHEDRON, allCells);


    /*
        ADD SCALAR DATA
    */
    // Define names for all scalar variables
    std::vector<std::string> scalarFieldNames = {"Pressure", "Velocity_X", "Velocity_Y", "Velocity_Z"};

    // Vector of vectors to store all scalar data for all 4 variables -> pressure, velX, velY, velZ
    std::vector<std::vector<double>> scalarValues(4, std::vector<double> (matrixSize));

    // Store data in scalarValues
    for (int i = 0; i < matrixSize; i++)
    {
        scalarValues[0][i] = pressure.variable(i) / listOfAllCells[i].getCellVolume(); // For pressure, divide by volume for output
        scalarValues[1][i] = velocity.variableX(i);
        scalarValues[2][i] = velocity.variableY(i);
        scalarValues[3][i] = velocity.variableZ(i);
    }

    // For every scalar variable, set cell centre values 
    for (int var = 0; var < 4; var++)
    {
        // Create a VTK data array for the scalar variable
        vtkSmartPointer<vtkDoubleArray> scalarData = vtkSmartPointer<vtkDoubleArray>::New();

        // Set the name of the variable
        scalarData->SetName(scalarFieldNames[var].c_str()); 

        // Set values for each cell
        for (int i = 0; i < matrixSize; i++)
        {
            scalarData->InsertNextValue(scalarValues[var][i]);
        }

        grid->GetCellData()->AddArray(scalarData);
    }


    /*
        ADD VECTOR DATA
    */
    vtkSmartPointer<vtkDoubleArray> velocityData = vtkSmartPointer<vtkDoubleArray>::New();
    velocityData->SetName("Velocity");
    velocityData->SetNumberOfComponents(3); // Assumes 3 vector components

    for (int i = 0; i < velocity.getMatrixSize(); i++)
    {
        std::array<double, 3> velData = {velocity.variableX(i), velocity.variableY(i), velocity.variableZ(i)};
        velocityData->InsertNextTuple(velData.data());
    }

    grid->GetCellData()->SetVectors(velocityData);


    /*
        Output VTK XML file
    */
    // Write the grid to a VTK XML unstructured grid file
    vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();

    std::stringstream ss;
    ss << filename + "SRF_Itr_" << iteration << ".vtu";
    std::string fileName = ss.str();
    writer->SetFileName(fileName.c_str()); 

    writer->SetDataModeToBinary();
    writer->SetInputData(grid);
    writer->Write(); 
}