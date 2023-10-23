CMAKE_MINIMUM_REQUIRED(VERSION 3.13.5)

project(TDSESolver)

find_package(OpenMP REQUIRED)

set(tdsesolver_SOURCES
    # C++
    src/cwrapper.cpp
    src/diagnostics.cpp
    src/fields.cpp
    src/hamiltonian.cpp
    src/hamiltonian_X.cpp
    src/hamiltonian_XZ.cpp
    src/hamiltonian_RZ.cpp
    src/hamiltonian_XYZ.cpp
    src/parameters.cpp
    src/potential.cpp
    src/probe.cpp
    src/probe_factory.cpp
    src/probe.X.cpp
    src/probe.XZ.cpp
    src/probe.RZ.cpp
    src/probe.XYZ.cpp
    src/tdsesolver.cpp
    src/tdsesolver_X.cpp
    src/tdsesolver_XZ.cpp
    src/tdsesolver_RZ.cpp
    src/tdsesolver_XYZ.cpp
    src/utils.cpp
    src/wavefunction.cpp
    src/wavefunction_X.cpp
    src/wavefunction_XZ.cpp
    src/wavefunction_RZ.cpp
    src/wavefunction_XYZ.cpp
)

if (TARGET propag-plugins)
    list(APPEND tdsesolver_SOURCES src/harmonic_TDSESolver.cpp)
    add_library(tdsesolver ${tdsesolver_SOURCES})
    target_link_libraries(tdsesolver PUBLIC propag)
    add_dependencies(tdsesolver propag)
else()
    add_library(tdsesolver SHARED ${tdsesolver_SOURCES})
endif()
target_include_directories(tdsesolver PUBLIC ${CMAKE_CURRENT_BINARY_DIR} include)
target_compile_features(tdsesolver PRIVATE cxx_std_17)
target_link_libraries(tdsesolver PUBLIC OpenMP::OpenMP_CXX)
set(CMAKE_CXX_FLAGS "-O3 -Wall -Wextra -ffast-math -g")
#target_link_libraries(tdsesolver PUBLIC propag)
file(COPY "${CMAKE_SOURCE_DIR}/pytdsesolver" DESTINATION "${CMAKE_BINARY_DIR}")
file(COPY "${CMAKE_SOURCE_DIR}/launch.py" DESTINATION "${CMAKE_BINARY_DIR}")
file(COPY "${CMAKE_SOURCE_DIR}/scan.py" DESTINATION "${CMAKE_BINARY_DIR}")
file(COPY "${CMAKE_SOURCE_DIR}/run.sh" DESTINATION "${CMAKE_BINARY_DIR}")

add_executable(tdsesolverapp src/main.cpp)
set_target_properties(tdsesolverapp PROPERTIES LINKER_LANGUAGE CXX)
target_link_libraries(tdsesolverapp PUBLIC tdsesolver)

