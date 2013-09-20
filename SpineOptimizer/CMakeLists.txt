cmake_minimum_required(VERSION 2.6)

PROJECT(spineOptimizer CXX)

# find_path( EIGEN_DIR Eigen DOC "A path to Eigen libaray, a linear algebra solver" )
# INCLUDE_DIRECTORIES( ${EIGEN_DIR} )

FIND_PACKAGE(ITK REQUIRED)
INCLUDE(${ITK_USE_FILE})

# FIND_PACKAGE(VTK REQUIRED)
# INCLUDE(${VTK_USE_FILE})

FIND_PACKAGE(Qt4)
INCLUDE(${QT_USE_FILE})
INCLUDE_DIRECTORIES(${QT_INCLUDE_DIR})

find_path( nlopt_DIR nlopt-2.3 DOC "A path to non-linear optimization library")
INCLUDE_DIRECTORIES(${nlopt_DIR}/api)

ADD_EXECUTABLE(spineOptimizer spineOptimizer.cpp)
# TARGET_LINK_LIBRARIES(spineOptimizer ${VTK_LIBRARIES})
TARGET_LINK_LIBRARIES(spineOptimizer ${ITK_LIBRARIES})
TARGET_LINK_LIBRARIES(spineOptimizer ${QT_LIBRARIES})
TARGET_LINK_LIBRARIES(spineOptimizer ${nlopt_DIR}/${CMAKE_CFG_INTDIR}/nlopt.lib)
