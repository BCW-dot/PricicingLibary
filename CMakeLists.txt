cmake_minimum_required(VERSION 3.12)
project(PricingLibrary VERSION 1.0.0 LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Build type
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# Compiler-specific flags for financial computing
if(MSVC)
    # MSVC compiler flags
    set(CMAKE_CXX_FLAGS_RELEASE "/O2 /Ob2 /DNDEBUG /fp:fast")
    set(CMAKE_CXX_FLAGS_DEBUG "/Od /Zi /DEBUG /fp:precise")
    add_compile_options(/W4 /permissive-)
else()
    # GCC/Clang compiler flags for high-performance financial computing
    set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -ffast-math -march=native -mtune=native")
    set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g -Wall -Wextra -Wpedantic")
    
    # Additional optimization flags for financial computing
    add_compile_options(
        -fno-signed-zeros          # Assume no signed zeros for faster floating-point
        -ffinite-math-only         # Assume finite math for optimization
        -funroll-loops             # Unroll loops for performance
        -ftree-vectorize           # Enable auto-vectorization
    )
endif()

# Find OpenMP for parallel computing (common in financial applications)
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    message(STATUS "OpenMP found - enabling parallel computing support")
else()
    message(STATUS "OpenMP not found - continuing without parallel support")
endif()

# Collect all source files by category
file(GLOB_RECURSE MODELS_SOURCES "src/models/*.cpp")
file(GLOB_RECURSE OPTIONS_SOURCES "src/options/*.cpp")
file(GLOB_RECURSE PRICING_SOURCES "src/pricing/*.cpp")
file(GLOB_RECURSE PORTFOLIO_SOURCES "src/portfolio/*.cpp")
file(GLOB_RECURSE UTILS_SOURCES "src/utils/*.cpp")
file(GLOB_RECURSE CHARTS_SOURCES "src/charts/*.cpp")
file(GLOB_RECURSE PDE_SOURCES "src/pde/*.cpp")

# Core library source files
set(PRICING_LIB_SOURCES
    src/Priceable.cpp
    src/testing.cpp
    ${MODELS_SOURCES}
    ${OPTIONS_SOURCES}
    ${PRICING_SOURCES}
    ${PORTFOLIO_SOURCES}
    ${UTILS_SOURCES}
    ${CHARTS_SOURCES}
    ${PDE_SOURCES}
)

# Collect all header files
file(GLOB_RECURSE MODELS_HEADERS "include/PricingLibrary/models/*.h")
file(GLOB_RECURSE OPTIONS_HEADERS "include/PricingLibrary/options/*.h")
file(GLOB_RECURSE PRICING_HEADERS "include/PricingLibrary/pricing/*.h")
file(GLOB_RECURSE PORTFOLIO_HEADERS "include/PricingLibrary/portfolio/*.h")
file(GLOB_RECURSE UTILS_HEADERS "include/PricingLibrary/utils/*.h")
file(GLOB_RECURSE CHARTS_HEADERS "include/PricingLibrary/charts/*.h")
file(GLOB_RECURSE PDE_HEADERS "include/PricingLibrary/pde/*.h")

set(PRICING_LIB_HEADERS
    include/PricingLibrary/Priceable.h
    include/PricingLibrary/stdafx.h
    include/PricingLibrary/testing.h
    ${MODELS_HEADERS}
    ${OPTIONS_HEADERS}
    ${PRICING_HEADERS}
    ${PORTFOLIO_HEADERS}
    ${UTILS_HEADERS}
    ${CHARTS_HEADERS}
    ${PDE_HEADERS}
)

# Create the pricing library
add_library(PricingLibrary STATIC ${PRICING_LIB_SOURCES} ${PRICING_LIB_HEADERS})

# Set library properties
set_target_properties(PricingLibrary PROPERTIES
    VERSION ${PROJECT_VERSION}
    SOVERSION 1
    ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib
)

# Include directories
target_include_directories(PricingLibrary
    PUBLIC 
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/include/PricingLibrary
)

# Link OpenMP if available
if(OpenMP_CXX_FOUND)
    target_link_libraries(PricingLibrary PUBLIC OpenMP::OpenMP_CXX)
endif()

# Precompiled headers for faster builds
target_precompile_headers(PricingLibrary PRIVATE include/PricingLibrary/stdafx.h)

# Demo executable
add_executable(PricingDemo examples/demo.cpp)
target_link_libraries(PricingDemo PRIVATE PricingLibrary)
set_target_properties(PricingDemo PROPERTIES
    OUTPUT_NAME "pricing_demo"
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
)

# Testing executable
add_executable(PricingTests EXCLUDE_FROM_ALL src/testing.cpp)
target_link_libraries(PricingTests PRIVATE PricingLibrary)
set_target_properties(PricingTests PROPERTIES
    OUTPUT_NAME "pricing_tests"
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
)

# Installation rules
include(GNUInstallDirs)

# Install library
install(TARGETS PricingLibrary
    EXPORT PricingLibraryTargets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/PricingLibrary
)

# Install headers preserving directory structure
install(DIRECTORY include/PricingLibrary/
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/PricingLibrary
    FILES_MATCHING PATTERN "*.h"
)

# Install demo executable
install(TARGETS PricingDemo
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

# Install examples
install(DIRECTORY examples/
    DESTINATION ${CMAKE_INSTALL_DOCDIR}/examples
    FILES_MATCHING PATTERN "*.cpp"
)

# Install CMake config files for find_package() support
install(EXPORT PricingLibraryTargets
    FILE PricingLibraryTargets.cmake
    NAMESPACE PricingLibrary::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/PricingLibrary
)

# Create config template if it doesn't exist
set(CONFIG_TEMPLATE_CONTENT "
include(CMakeFindDependencyMacro)

find_dependency(OpenMP)

include(\"\${CMAKE_CURRENT_LIST_DIR}/PricingLibraryTargets.cmake\")

check_required_components(PricingLibrary)
")

file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/PricingLibraryConfig.cmake.in" "${CONFIG_TEMPLATE_CONTENT}")

# Create config file
include(CMakePackageConfigHelpers)
configure_package_config_file(
    "${CMAKE_CURRENT_BINARY_DIR}/PricingLibraryConfig.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/PricingLibraryConfig.cmake"
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/PricingLibrary
)

write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/PricingLibraryConfigVersion.cmake"
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion
)

install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/PricingLibraryConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/PricingLibraryConfigVersion.cmake"
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/PricingLibrary
)

# Custom targets
add_custom_target(run
    COMMAND $<TARGET_FILE:PricingDemo>
    DEPENDS PricingDemo
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/bin
    COMMENT "Running pricing library demo"
)

add_custom_target(tests
    COMMAND $<TARGET_FILE:PricingTests>
    DEPENDS PricingTests
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/bin
    COMMENT "Running pricing library tests"
)

# Create a comprehensive summary target
add_custom_target(info
    COMMAND ${CMAKE_COMMAND} -E echo "PricingLibrary v${PROJECT_VERSION} - Quantitative Finance Library"
    COMMAND ${CMAKE_COMMAND} -E echo ""
    COMMAND ${CMAKE_COMMAND} -E echo "Directory Structure:"
    COMMAND ${CMAKE_COMMAND} -E echo "  include/PricingLibrary/  - Public headers"
    COMMAND ${CMAKE_COMMAND} -E echo "    models/               - Financial models (Black-Scholes, etc.)"
    COMMAND ${CMAKE_COMMAND} -E echo "    options/              - Options classes (Call, Put, Asian, etc.)"
    COMMAND ${CMAKE_COMMAND} -E echo "    pricing/              - Pricing engines (Monte Carlo, etc.)"
    COMMAND ${CMAKE_COMMAND} -E echo "    portfolio/            - Portfolio management & hedging"
    COMMAND ${CMAKE_COMMAND} -E echo "    utils/                - Mathematical utilities"
    COMMAND ${CMAKE_COMMAND} -E echo "    pde/                  - PDE and numerical methods"
    COMMAND ${CMAKE_COMMAND} -E echo "    charts/               - Visualization tools"
    COMMAND ${CMAKE_COMMAND} -E echo ""
    COMMAND ${CMAKE_COMMAND} -E echo "Build Targets:"
    COMMAND ${CMAKE_COMMAND} -E echo "  PricingLibrary - Static library"
    COMMAND ${CMAKE_COMMAND} -E echo "  PricingDemo    - Demo executable"
    COMMAND ${CMAKE_COMMAND} -E echo "  PricingTests   - Test executable"
    COMMAND ${CMAKE_COMMAND} -E echo ""
    COMMAND ${CMAKE_COMMAND} -E echo "Commands:"
    COMMAND ${CMAKE_COMMAND} -E echo "  make run       - Run demo"
    COMMAND ${CMAKE_COMMAND} -E echo "  make tests     - Run tests"
    COMMAND ${CMAKE_COMMAND} -E echo "  make install   - Install library"
    COMMENT "Displaying PricingLibrary information"
)

# Print build information
message(STATUS "")
message(STATUS "PricingLibrary v${PROJECT_VERSION} Configuration:")
message(STATUS "  Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "  C++ standard: ${CMAKE_CXX_STANDARD}")
message(STATUS "  Compiler: ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")
message(STATUS "  OpenMP support: ${OpenMP_CXX_FOUND}")
message(STATUS "")
message(STATUS "Project Structure:")
message(STATUS "  Source files: ${PRICING_LIB_SOURCES}")
message(STATUS "  Header files: ${PRICING_LIB_HEADERS}")
message(STATUS "")
message(STATUS "Available targets: PricingLibrary, PricingDemo, PricingTests, run, tests, info")
message(STATUS "")