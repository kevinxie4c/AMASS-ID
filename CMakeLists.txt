project(inverse)

find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIRS})

file(GLOB EXT_SRC_FILES ${PROJECT_SOURCE_DIR}/src/ext/*.cpp)
file(GLOB BUTTER_SRC_FILES ${PROJECT_SOURCE_DIR}/src/c_butterworth/*.cpp)
set (CMAKE_CXX_FLAGS "-std=c++14")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS} -g")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} -O3 -DNDEBUG")
add_executable(inverse.out src/inverse.cpp src/BVHData.cpp ${EXT_SRC_FILES} ${BUTTER_SRC_FILES})
target_link_libraries(inverse.out m stdc++ dart dart-collision-ode assimp boost_system)
