PROJECT(detector)

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --std=c++98 -W -Wall -Wno-sign-compare -O3 -s -pipe -Isrc")

SET(SOURCES
  src/main.cpp
)

ADD_EXECUTABLE(detector ${SOURCES})
TARGET_LINK_LIBRARIES(detector)
