cmake_minimum_required(VERSION 3.9)

# define project name, version
project(PSEMolDyn_Group1 VERSION 0.0.1)

#option to enable compilation with linter
option(SANITIZE "Enable compilation with sanitizer" OFF)

option(TESTS "Enable generation of tests" ON)

option(PROFILE "Enable compilation with symbols" OFF)

option(VECTOR "Enable vectorization" OFF)

option(FAST_I "Enable fastet Benchmark" OFF)

option(OMP "Enable OpenMP" OFF)

option(IOMP "Enable OpenMP for icpc" OFF)


option(FAST_G "Enable fastet Benchmark GCC" OFF)
# option to disable logging in performance critical regions
option(BMARK "Disable io in benchmark region" OFF)
# let ccmake and cmake-gui offer the default build type options
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug;Release;RelWithDebInfo;MinSizeRel")

# set Release as the default build type if it is not yet set.
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING
            "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif (NOT CMAKE_BUILD_TYPE)

# collect all cpp files
file(GLOB_RECURSE MY_SRC
        "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp"
        # header don't need to be included but this might be necessary for some IDEs
        "${CMAKE_CURRENT_SOURCE_DIR}/src/*.h"
        )
if (FAST_I)
    message("Fast ICPC")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -ipo -fast")
    SET(CMAKE_RANLIB ${CMAKE_CXX_COMPILER_RANLIB} CACHE PATH "RANLIB" FORCE)
endif ()

if (OMP)
    message("with openmp")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
endif ()

if (IOMP)
    message("with openmp")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -qopenmp")
endif ()

if (FAST_G)
    message("FAST G++")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Ofast")
    SET(CMAKE_RANLIB ${CMAKE_CXX_COMPILER_RANLIB} CACHE PATH "RANLIB" FORCE)
endif ()

# create make target
add_executable(MolSim ${MY_SRC} src/Statistics.cpp src/Statistics.h)
# set cxx standard. You may raise this if you want.
target_compile_features(MolSim
        PRIVATE
        cxx_std_17
        )

if (PROFILE)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg -O2")
    SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg")
    SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -pg")
endif ()

if (VECTOR)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
    target_compile_definitions(MolSim PRIVATE AVX)
endif ()

if (BMARK)
    target_compile_definitions(MolSim PRIVATE SPDLOG_ACTIVE_LEVEL=SPDLOG_LEVEL_OFF BENCHMARK)
endif ()
target_include_directories(MolSim
        PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/libs/libxsd
        PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src
        )

find_package(PkgConfig REQUIRED)
pkg_check_modules(XERCESC REQUIRED IMPORTED_TARGET xerces-c)
target_link_libraries(MolSim
        # stuff that is used in headers and source files
        PUBLIC
        PkgConfig::XERCESC
        PUBLIC
        spdlog
        )

# if the option is set, the program will be compiled with fsanitize=address
if (SANITIZE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address")
endif ()

# activate all compiler warnings. Clean up your code :P
# depending on the compiler different flags are used
target_compile_options(MolSim
        PRIVATE
        $<$<CXX_COMPILER_ID:GNU>:-Wall>
        # disable some spam warnings for icpc...
        $<$<CXX_COMPILER_ID:Intel>:-w3 -wd383,981,1418,1572,2259>
        )


list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/modules)

include(doxygen)
include(spdlog)

if (TESTS)
    #set up google test
    include(FetchContent)

    FetchContent_Declare(
            googletest
            GIT_REPOSITORY https://github.com/google/googletest.git
            GIT_TAG 58d77fa8070e8cec2dc1ed015d66b454c8d78850
    )

    FetchContent_MakeAvailable(googletest)

    enable_testing()

    file(GLOB_RECURSE MY_TEST
            "${CMAKE_CURRENT_SOURCE_DIR}/tests/*.cc"
            )
    add_executable(
            MolSimTest
            ${MY_TEST}
    )
    file(GLOB_RECURSE FORCE
            "${CMAKE_CURRENT_SOURCE_DIR}/src/forceCalculation/*.cpp"
            )
    file(GLOB_RECURSE INPUT
            "${CMAKE_CURRENT_SOURCE_DIR}/src/inputReader/*.cpp"
            "${CMAKE_CURRENT_SOURCE_DIR}/src/utils/*.cpp"
            )

    file(GLOB_RECURSE CONTAINER
            "${CMAKE_CURRENT_SOURCE_DIR}/src/container/*.cpp"
            )


    add_library(lib_test
            ${FORCE}
            ${INPUT}
            ${CMAKE_CURRENT_SOURCE_DIR}/src/utils/ArrayUtils.h
            ${CMAKE_CURRENT_SOURCE_DIR}/src/MolSimLogger.cpp
            ${CONTAINER}
            ${CMAKE_CURRENT_SOURCE_DIR}/src/Reflecting.cpp
            ${CMAKE_CURRENT_SOURCE_DIR}/src/Particle.cpp
            ${CMAKE_CURRENT_SOURCE_DIR}/src/outputWriter/FileWriter.cpp
            ${CMAKE_CURRENT_SOURCE_DIR}/src/Simulation.cpp
            ${CMAKE_CURRENT_SOURCE_DIR}/src/utils/Thermostat.cpp
            ${CMAKE_CURRENT_SOURCE_DIR}/src/Statistics.cpp src/utils/BoundaryException.cpp src/utils/BoundaryException.h)


    target_link_libraries(lib_test PRIVATE spdlog::spdlog)
    target_compile_features(lib_test
            PRIVATE
            cxx_std_20
            )
    target_include_directories(lib_test
            PRIVATE
            ${CMAKE_CURRENT_SOURCE_DIR}/src/
            PUBLIC
            ${CMAKE_CURRENT_SOURCE_DIR}/libs/libxsd/)
    target_include_directories(MolSimTest
            PRIVATE
            ${CMAKE_CURRENT_SOURCE_DIR}/src/
            )
    target_link_libraries(MolSimTest
            GTest::gtest
            gmock
            lib_test
            PkgConfig::XERCESC
            spdlog::spdlog
            )

    include(GoogleTest)

    gtest_discover_tests(MolSimTest)
endif ()






