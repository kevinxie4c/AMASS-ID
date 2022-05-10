project(inverse)
cmake_minimum_required(VERSION 3.16)

find_package(Eigen3 REQUIRED)
#find_package(Python3 REQUIRED)	# Uncomment this if you want to use system python

include_directories(${EIGEN3_INCLUDE_DIRS})

file(GLOB EXT_SRC_FILES ${PROJECT_SOURCE_DIR}/src/ext/*.cpp)
file(GLOB BUTTER_SRC_FILES ${PROJECT_SOURCE_DIR}/src/c_butterworth/*.cpp)

execute_process(COMMAND python3-config --prefix OUTPUT_VARIABLE Python3_PREFIX OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND python3-config --includes OUTPUT_VARIABLE Python3_INCLUDE_DIRS OUTPUT_STRIP_TRAILING_WHITESPACE)
#execute_process(COMMAND python3-config --ldflags OUTPUT_VARIABLE Python3_LINK_OPTIONS OUTPUT_STRIP_TRAILING_WHITESPACE)

separate_arguments(Python3_INCLUDE_DIRS UNIX_COMMAND "${Python3_INCLUDE_DIRS}")
#separate_arguments(Python3_LINK_OPTIONS UNIX_COMMAND "${Python3_LINK_OPTIONS}")
add_compile_options(${Python3_INCLUDE_DIRS})

set(Python3_LINK_OPTIONS -L/home/kevin/miniconda3/envs/py3.7/lib -Wl,-rpath=${Python3_PREFIX}/lib -lpython3.7m -lcrypt -ldl -lutil -lm) # Use python in the current environment. Comment this if you want to use system python

message(STATUS "Python3_PREFIX: ${Python3_PREFIX}")
message(STATUS "Python3_INCLUDE_DIRS: ${Python3_INCLUDE_DIRS}")
message(STATUS "Python3_LINK_OPTIONS: ${Python3_LINK_OPTIONS}")

set(CMAKE_CXX_FLAGS "-std=c++14")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS} -g")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} -O3 -DNDEBUG")

add_executable(inverse.out src/inverse.cpp src/SimCharacter.cpp src/PyUtil.cpp src/IOUtil.cpp ${EXT_SRC_FILES} ${BUTTER_SRC_FILES})
target_link_libraries(inverse.out m stdc++ dart dart-collision-ode assimp boost_system ${Python3_LINK_OPTIONS})

add_executable(CharacterViewer.out src/CharacterViewer.cpp src/SimCharacter.cpp src/PyUtil.cpp src/IOUtil.cpp)
target_link_libraries(CharacterViewer.out stdc++ dart dart-gui-osg osg osgViewer osgUtil assimp boost_system dart-external-imgui ${Python3_LINK_OPTIONS})
