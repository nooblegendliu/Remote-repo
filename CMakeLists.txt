cmake_minimum_required(VERSION 3.25)
project(Repo VERSION 0.1.0 LANGUAGES C CXX)

include(CTest)
enable_testing()

file(GLOB SOURCES "src/*.cpp")
add_executable(Repo ${SOURCES})
target_include_directories(Repo PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/inc)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
