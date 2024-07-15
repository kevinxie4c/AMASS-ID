project(inverse)
cmake_minimum_required(VERSION 3.16)

find_package(Eigen3 REQUIRED)
find_package(Python3 REQUIRED)	# Uncomment this if you want to use system python
find_package(DART REQUIRED COMPONENTS collision-ode gui-osg)
find_package(nlohmann_json REQUIRED)
find_package(OpenSceneGraph REQUIRED)

include_directories(${EIGEN3_INCLUDE_DIRS})

file(GLOB EXT_SRC_FILES ${PROJECT_SOURCE_DIR}/src/ext/*.cpp)
file(GLOB BUTTER_SRC_FILES ${PROJECT_SOURCE_DIR}/src/c_butterworth/*.cpp)

execute_process(COMMAND python3-config --prefix OUTPUT_VARIABLE Python3_PREFIX OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND python3-config --includes OUTPUT_VARIABLE Python3_INCLUDE_DIRS OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND python3-config --ldflags OUTPUT_VARIABLE Python3_LINK_OPTIONS OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND perl -E "`python3 --version`=~/(3\\.\\d*)/; say $1" OUTPUT_VARIABLE Python3_VERSION OUTPUT_STRIP_TRAILING_WHITESPACE)

separate_arguments(Python3_INCLUDE_DIRS UNIX_COMMAND "${Python3_INCLUDE_DIRS}")
separate_arguments(Python3_LINK_OPTIONS UNIX_COMMAND "${Python3_LINK_OPTIONS}")
add_compile_options(${Python3_INCLUDE_DIRS})

message(STATUS "Python3_PREFIX: ${Python3_PREFIX}")
message(STATUS "Python3_INCLUDE_DIRS: ${Python3_INCLUDE_DIRS}")
message(STATUS "Python3_LINK_OPTIONS: ${Python3_LINK_OPTIONS}")
message(STATUS "Python3_VERSION: ${Python3_VERSION}")

# Sometimes python3-config output flags that miss some pieces, then you can uncomment the next line.
# set(PYTHON_LD_FLAGS_EXTRA "python${Python3_VERSION} -Wl,-rpath=${Python3_PREFIX}/lib")

set(CMAKE_CXX_FLAGS "-std=c++14")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS} -g")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} -O3 -DNDEBUG")

if(DEFINED MOSEK_INCLUDE_PATH)
    include_directories(${MOSEK_INCLUDE_PATH})
else()
    message(WARNING "You might need to define MOSEK_INCLUDE_PATH")
endif()

if(DEFINED MOSEK_LIB_PATH)
    set(MOSEK_LINK_OPTIONS -L${MOSEK_LIB_PATH} -Wl,-rpath=${MOSEK_LIB_PATH} -lmosek64)
    message(STATUS "${MOSEK_LINK_OPTIONS}")
else()
    message(WARNING "You might need to define MOSEK_LIB_PATH")
endif()


add_executable(inverse.out src/inverse.cpp src/SimCharacter.cpp src/PyUtil.cpp src/IOUtil.cpp ${EXT_SRC_FILES} ${BUTTER_SRC_FILES})
target_link_libraries(inverse.out m stdc++ dart dart-collision-ode assimp boost_system ${Python3_LINK_OPTIONS} ${PYTHON_LD_FLAGS_EXTRA} ${MOSEK_LINK_OPTIONS})

add_executable(CharacterViewer.out src/CharacterViewer.cpp src/SimCharacter.cpp src/PyUtil.cpp src/IOUtil.cpp)
target_link_libraries(CharacterViewer.out stdc++ dart dart-gui-osg osg osgViewer osgUtil osgGA assimp boost_system dart-external-imgui ${PYTHON_LD_FLAGS_EXTRA} ${Python3_LINK_OPTIONS})
