cmake_minimum_required(VERSION 2.6)

PROJECT(QuadEdge CXX)

find_path( EIGEN_DIR Eigen DOC "A path to Eigen libaray, a linear algebra solver" )
INCLUDE_DIRECTORIES( ${EIGEN_DIR} )

SET(sources cell.cc edge.cc face.cc obj.cc vertex.cc)
SET(headers array.hh cell.hh edge.hh face.hh list.hh obj.hh vertex.hh)
 
ADD_LIBRARY(QuadEdge STATIC ${sources} ${headers})
