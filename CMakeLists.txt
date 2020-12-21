cmake_minimum_required(VERSION 3.10)
project(EFGE)

set(CMAKE_CXX_STANDARD 14)

include_directories(lib)

set(SOURCE_FILES src/Graph.cpp src/Model.cpp src/Unigram.cpp src/Vocabulary.cpp src/Node.cpp src/Utilities.cpp)
set(LIB_FILES lib/Graph.h lib/Model.h lib/Unigram.h lib/Vocabulary.h lib/Node.h lib/Utilities.h)

add_executable(${CMAKE_PROJECT_NAME} main.cpp ${SOURCE_FILES} ${LIB_FILES})
