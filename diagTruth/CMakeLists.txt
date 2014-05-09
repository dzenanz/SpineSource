cmake_minimum_required(VERSION 2.6)

PROJECT(diagTruth CXX)

INCLUDE_DIRECTORIES( ${CMAKE_CURRENT_SOURCE_DIR} )
add_subdirectory(qe)
INCLUDE_DIRECTORIES(qe)

find_path(lp_solve_DIR lp_solve_5.5 DOC "A path to lpsolve.sourceforge.net/5.5/")
INCLUDE_DIRECTORIES( ${lp_solve_DIR} )
find_library(lp_solve_LIB lp_solve_5.5 DOC "A path to debug/release containing folder")

FIND_PACKAGE(ITK REQUIRED)
INCLUDE(${ITK_USE_FILE})

FIND_PACKAGE(VTK REQUIRED)
INCLUDE(${VTK_USE_FILE})

find_path( EIGEN_DIR Eigen DOC "A path to Eigen libaray, a linear algebra solver" )
INCLUDE_DIRECTORIES( ${EIGEN_DIR} )

ADD_EXECUTABLE(diagTruth diagTruth.cpp)

TARGET_LINK_LIBRARIES(diagTruth ${VTK_LIBRARIES})
TARGET_LINK_LIBRARIES(diagTruth ${ITK_LIBRARIES})
TARGET_LINK_LIBRARIES(diagTruth ${lp_solve_LIB}/${CMAKE_CFG_INTDIR}/liblpsolve55.lib)
TARGET_LINK_LIBRARIES(diagTruth QuadEdge)