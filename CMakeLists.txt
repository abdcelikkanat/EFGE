cmake_minimum_required(VERSION 3.12)
project(fast_bern)

set(CMAKE_CXX_STANDARD 14)

include_directories(.)

add_executable(fast_bern
        Graph.cpp
        Graph.h
        main.cpp
        Model.cpp
        Model.h
        Unigram.cpp
        Unigram.h Vocabulary.cpp Vocabulary.h Node.cpp Node.h Model.cpp Model.h)
