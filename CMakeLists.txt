cmake_minimum_required(VERSION 3.16)
project(example)

list(PREPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

add_definitions(-D _USE_MATH_DEFINES)

# Libigl
option(LIBIGL_WITH_OPENGL            "Use OpenGL"         ON)
option(LIBIGL_WITH_OPENGL_GLFW       "Use GLFW"           ON)
include(libigl)
include_directories(D:/OpenMesh/OpenMesh/include)

link_directories(D:/OpenMesh/OpenMesh/lib)

add_subdirectory(src)


