cmake_minimum_required(VERSION 2.8.12.2)

INCLUDE_DIRECTORIES( ${CMAKE_CURRENT_SOURCE_DIR} )
add_subdirectory(qe)
INCLUDE_DIRECTORIES(qe)

PROJECT(Spine CXX C)

SET(CMAKE_PREFIX_PATH "C:/Program Files (x86)/Windows Kits/8.0/Lib/win8/um/x64;C:/Program Files (x86)/Microsoft SDKs/Windows/v7.0A/Lib")

find_path( EIGEN_DIR Eigen DOC "A path to Eigen libaray, eigen.tuxfamily.org" )
INCLUDE_DIRECTORIES( ${EIGEN_DIR} )

#find_path( LevMar_DIR levmar-2.6 DOC "A path to Levenberg-Marquardt library")
#INCLUDE_DIRECTORIES( ${LevMar_DIR} )

find_path(lp_solve_DIR lp_solve_5.5 DOC "A path to lpsolve.sourceforge.net/5.5/")
INCLUDE_DIRECTORIES( ${lp_solve_DIR} )
find_library(lp_solve_LIB lp_solve_5.5 DOC "A path to debug/release containing folder") #D:/SDK/lp_solve_5.5/lpsolve55/x64/


FIND_PACKAGE(OpenCV REQUIRED)
#message("OpenCV_LIBS=${OpenCV_LIBS}")
FIND_PACKAGE(TBB REQUIRED) #required by OpenCV for threaded object recognition
include_directories(${TBB_INCLUDE_DIRS})
link_directories(${TBB_LIBRARY_DIRS})

FIND_PACKAGE(ITK REQUIRED)
INCLUDE(${ITK_USE_FILE})

FIND_PACKAGE(VTK REQUIRED COMPONENTS vtkRenderingVolumeOpenGL vtkGUISupportQt vtkGUISupportQtOpenGL vtkRenderingQt vtkRenderingOpenGL NO_MODULE)
#FIND_PACKAGE(VTK REQUIRED)
INCLUDE(${VTK_USE_FILE})

#set(CMAKE_AUTOMOC TRUE)
FIND_PACKAGE(Qt5Gui)
FIND_PACKAGE(Qt5OpenGL)
#FIND_PACKAGE(Qt5Concurrent)
   
SET(SAMPLE_SRCS
    Spine.cxx
    mainWindow.cpp
    slicePainter.cpp
    itkDatImageIO.cxx
    itkDatImageIO.h
    itkDatImageIOFactory.cxx
    itkDatImageIOFactory.h
    vec3.h
    MainLogic.cpp
    declarations.h
    butterfly.cpp
    qeUtil.cpp
    LHstuff.cpp
    itkStructureTensorImageFilter.h
    itkStructureTensorImageFilter.txx
    itkCustomCanny.h
    itkCustomCanny.txx
    voxelClassification.h
    voxelClassification.cpp
    vertebra.h
    vertebra.cpp
    polyfit.cpp
    vtkOBJWriter.cxx
    )
SET(SAMPLE_MOC_HDRS mainWindow.h MainLogic.h slicePainter.h)
SET(SAMPLE_UIS mainWindow.ui slicePainter.ui)
SET(SAMPLE_RCS mainWindow.qrc)

QT5_ADD_RESOURCES( SAMPLE_RC_SRCS ${SAMPLE_RCS} )
QT5_WRAP_UI( SAMPLE_UI_HDRS ${SAMPLE_UIS} )
QT5_WRAP_CPP( SAMPLE_MOC_SRCS ${SAMPLE_MOC_HDRS} )

ADD_EXECUTABLE(Spine ${SAMPLE_SRCS} ${SAMPLE_MOC_HDRS} ${SAMPLE_MOC_SRCS} ${SAMPLE_RC_SRCS} ${SAMPLE_UI_HDRS})

TARGET_LINK_LIBRARIES(Spine ${OpenCV_LIBS})
#TARGET_LINK_LIBRARIES(Spine ${LevMar_DIR}/${CMAKE_CFG_INTDIR}/levmar.lib)
TARGET_LINK_LIBRARIES(Spine ${lp_solve_LIB}/${CMAKE_CFG_INTDIR}/liblpsolve55.lib)
qt5_use_modules(Spine Gui OpenGL Concurrent)
TARGET_LINK_LIBRARIES(Spine ${ITK_LIBRARIES})
TARGET_LINK_LIBRARIES(Spine ${VTK_LIBRARIES})
TARGET_LINK_LIBRARIES(Spine QuadEdge)