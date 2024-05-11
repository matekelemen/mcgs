cmake_minimum_required(VERSION 3.15.0)
project(mcgs LANGUAGES CXX VERSION 0.0.1)

# Language requirements
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Options
string(TOUPPER ${PROJECT_NAME} PROJECT_NAME_UPPER)
option(${PROJECT_NAME_UPPER}_BUILD_TESTS "Compile tests." OFF)

# Find dependencies
find_package(OpenMP)

# Define library
add_library(${PROJECT_NAME} SHARED)
target_sources(${PROJECT_NAME} PRIVATE
               "${CMAKE_CURRENT_SOURCE_DIR}/src/color.cpp"
               "${CMAKE_CURRENT_SOURCE_DIR}/src/partition.cpp"
               "${CMAKE_CURRENT_SOURCE_DIR}/src/solve.cpp")
target_include_directories(${PROJECT_NAME}
                           PUBLIC "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
                           INTERFACE "$<INSTALL_INTERFACE:${PROJECT_NAME}/include>"
                                     "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>")
set_target_properties(${PROJECT_NAME} PROPERTIES
                      LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
                      ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
                      RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
                      INSTALL_RPATH "$ORIGIN")

if (OpenMP_CXX_FOUND)
    target_link_libraries(${PROJECT_NAME} PRIVATE OpenMP::OpenMP_CXX)
endif()

if (CMAKE_CXX_COMPILER_IF MATCHES "GNU|Clang")
    target_compile_options(${PROJECT_NAME} PRIVATE
                           -Wall -Wpedantic -Wextra -Werror)
endif()

# Define testrunner
if (${${PROJECT_NAME_UPPER}_BUILD_TESTS})
    add_executable(${PROJECT_NAME}_testrunner
                   "${CMAKE_CURRENT_SOURCE_DIR}/test/main.cpp"
                   "${CMAKE_CURRENT_SOURCE_DIR}/test/parseMatrixMarket.cpp")
    target_link_libraries(${PROJECT_NAME}_testrunner PRIVATE ${PROJECT_NAME})
    set_target_properties(${PROJECT_NAME}_testrunner PROPERTIES
                          RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
endif()

# Package
include(InstallRequiredSystemLibraries)
set(CPACK_PACKAGE_VENDOR "Máté Kelemen")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "A simple implementation of multicolor Gauss-Seidel relaxation.")
set(CPACK_PACKAGE_VERSION_MAJOR ${PROJECT_VERSION_MAJOR})
set(CPACK_PACKAGE_VERSION_MINOR ${PROJECT_VERSION_MINOR})
set(CPACK_PACKAGE_VERSION_PATCH ${PROJECT_VERSION_PATCH})
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/license")
set(CPACK_RESOURCE_FILE_README "${CMAKE_CURRENT_SOURCE_DIR}/readme")

set(CPACK_GENERATOR "ZIP")
set(CPACK_SOURCE_GENERATOR "ZIP")
set(CPACK_SOURCE_IGNORE_FILES
    /.git
    /build
    /.gitignore)
include(CPack)

# Generate CMake config
include(CMakePackageConfigHelpers)
configure_package_config_file("${CMAKE_CURRENT_SOURCE_DIR}/Config.cmake.in"
                              "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
                              INSTALL_DESTINATION "lib/cmake/${PROJECT_NAME}"
                              NO_SET_AND_CHECK_MACRO
                              NO_CHECK_REQUIRED_COMPONENTS_MACRO)
write_basic_package_version_file("${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
                                 VERSION "${${PROJECT_NAME}_VERSION_MAJOR}.${${PROJECT_NAME}_VERSION_MINOR}"
                                 COMPATIBILITY AnyNewerVersion)

# Install
install(TARGETS ${PROJECT_NAME}
        EXPORT ${PROJECT_NAME}Targets
        LIBRARY DESTINATION "lib"
        ARCHIVE DESTINATION "lib"
        RUNTIME DESTINATION "bin"
        INCLUDES DESTINATION "include")
install(EXPORT ${PROJECT_NAME}Targets
        FILE ${PROJECT_NAME}Targets.cmake
        NAMESPACE ${PROJECT_NAME}::
        DESTINATION "lib/cmake/${PROJECT_NAME}")
