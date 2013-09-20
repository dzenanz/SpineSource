cmake_minimum_required(VERSION 2.6)

PROJECT(centerRecTrainer CXX)

FIND_PACKAGE(ITK REQUIRED)
INCLUDE(${ITK_USE_FILE})

ADD_EXECUTABLE(centerRecTrainer main.cxx)

TARGET_LINK_LIBRARIES(centerRecTrainer ${ITK_LIBRARIES})
