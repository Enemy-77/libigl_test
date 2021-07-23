# Add your project files
set(SRC_FILES ${PROJECT_SOURCE_DIR}/src/main.cc ${PROJECT_SOURCE_DIR}/src/mesh.cc ${PROJECT_SOURCE_DIR}/src/mesh.h)
add_executable(${PROJECT_NAME} ${SRC_FILES})
target_link_libraries(${PROJECT_NAME} PUBLIC igl::core igl::opengl_glfw OpenMeshCored OpenMeshToolsd)
