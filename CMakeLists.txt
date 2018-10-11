cmake_minimum_required(VERSION 3.7.0)
project(dvr VERSION 0.1.0 LANGUAGES CXX)

file(GLOB sources ${CMAKE_CURRENT_SOURCE_DIR}/sources/*.cpp)
file(GLOB headers ${CMAKE_CURRENT_SOURCE_DIR}/sources/*.h)

add_executable(${PROJECT_NAME} ${sources} ${headers})

target_include_directories(${PROJECT_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/sources)

target_compile_features(${PROJECT_NAME} PRIVATE cxx_std_11)
target_compile_options(${PROJECT_NAME} PRIVATE -Wall -pipe)

find_package(Eigen3 3.3 REQUIRED NO_MODULE)
target_link_libraries(${PROJECT_NAME} Eigen3::Eigen)
