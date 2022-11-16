MolSim
===

# Group 1 #
Members:
* Lukas Neef
* Dominik Rammelsberger

# Code #
* Link: https://github.com/lneef/PSEMolDyn_Group1
* Branch: main
* Revision: 

# Run Instructions #
* Compiler: GNU 12.2.0, 11.3.0
* run mkdir build and cd build/
* run cmake ..
* run make or make all
* run ./MolSim --help for further information on usage or see Usage
* run ctest for running all tests

# Doxygen #
* Per default cmake .. creates a target for generating a doxygen documentation
* you can disable this behaviour by passing DISABLE\_DOXYGEN=ON
* to create the doxygen documentation run make doc\_doxygen

# CMake #
* by setting LINTER=ON compilation with fsanitize=address is enables
* by setting BMARK=ON compilation in benchmark mode is enabled, i.e. all io in benchmark region is disabled

# Usage #
MolSim is a simulator for simulating the movement of planets and molecules.

The following command line options are supported:

    -e      Expects an argument of type double. Using this option the end time of the simulation can be set.
            The default end time is 1000

    -f      Expects an argument of type double. Using this option the time step (delta\_t) of the simulation can be set.
            The default time step is 0.014

    -l     Expects as argument either off, info, debug, trace or error. Using this option the global log_level can be set

    --planet Usage: --planet path-to-file. Using this command line option particles are read from the specified file and
             the planet simulation is started.

    --cub   Using this option the molecular simulation is started. Particle are generated from cuboid which are either read
            from a file or from the command line. In the former case --cub=path-to-file has to be used and in the latter
            --cub.

    --help  This page is printed.