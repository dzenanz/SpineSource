cmake_minimum_required(VERSION 2.6)

PROJECT(gcCompare CXX)

FIND_PACKAGE(ITK REQUIRED)
INCLUDE(${ITK_USE_FILE})

FIND_PACKAGE(VTK REQUIRED)
INCLUDE(${VTK_USE_FILE})

ADD_EXECUTABLE(gcCompare main.cxx
    itkPGMImageIO.cxx itkPGMImageIO.h
    itkPGMImageIOFactory.cxx itkPGMImageIOFactory.h
    vtkOBJWriter.cxx vtkOBJWriter.h)

TARGET_LINK_LIBRARIES(gcCompare ${ITK_LIBRARIES})
TARGET_LINK_LIBRARIES(gcCompare ${VTK_LIBRARIES})
