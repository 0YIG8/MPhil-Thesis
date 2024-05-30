find_package(Eigen3 REQUIRED)
find_package(VTK REQUIRED)
find_package(PkgConfig REQUIRED)
pkg_search_module(LIBCONFIG++ REQUIRED libconfig++)

include_directories(/usr/include/vtk-9.1)

set(CMAKE_MODULE_PATH "usr/lib/x86_64-linux-gnu/cmake/vtk-9.1")


# For Lid-Driven Cavity
add_executable(main_Cavity Vector.C Point.C Face.C Cell.C BoundaryPatch.C Discretization.C Mesh.C MRFMesh.C SRFMesh.C ScalarField.C VectorField.C FileParser.C OutputData.C OutputVTK.C main_Cavity.C)
target_include_directories(main_Cavity PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_link_libraries(main_Cavity PUBLIC Eigen3::Eigen ${VTK_LIBRARIES} ${LIBCONFIG++_LINK_LIBRARIES})
target_include_directories(main_Cavity PUBLIC ${LIBCONFIG++_INCLUDE_DIRS})


# For SRF
add_executable(main_SRF Vector.C Point.C Face.C Cell.C BoundaryPatch.C MovingMeshDiscretization.C Mesh.C MRFMesh.C SRFMesh.C SIMPLE.C ScalarField.C VectorField.C FileParser.C OutputData.C OutputVTK.C main_SRF.C)
target_include_directories(main_SRF PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_link_libraries(main_SRF PUBLIC Eigen3::Eigen ${VTK_LIBRARIES} ${LIBCONFIG++_LINK_LIBRARIES})
target_include_directories(main_SRF PUBLIC ${LIBCONFIG++_INCLUDE_DIRS})


# For MRF
add_executable(main_MRF Vector.C Point.C Face.C Cell.C BoundaryPatch.C MovingMeshDiscretization.C Mesh.C MRFMesh.C SRFMesh.C SIMPLE.C ScalarField.C VectorField.C FileParser.C OutputData.C OutputVTK.C main_MRF.C)
target_include_directories(main_MRF PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_link_libraries(main_MRF PUBLIC Eigen3::Eigen ${VTK_LIBRARIES} ${LIBCONFIG++_LINK_LIBRARIES})
target_include_directories(main_MRF PUBLIC ${LIBCONFIG++_INCLUDE_DIRS})


# For test file
add_executable(main_test Vector.C Point.C Face.C Cell.C BoundaryPatch.C MovingMeshDiscretization.C Mesh.C MRFMesh.C SRFMesh.C SIMPLE.C ScalarField.C VectorField.C FileParser.C OutputData.C OutputVTK.C main_test.C)
target_include_directories(main_test PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_link_libraries(main_test PUBLIC Eigen3::Eigen ${VTK_LIBRARIES} ${LIBCONFIG++_LINK_LIBRARIES})
target_include_directories(main_test PUBLIC ${LIBCONFIG++_INCLUDE_DIRS})