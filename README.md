MolSim
===

# Group A #
Members:
* Lukas Neef
* Dominik Rammelsberger
* Xiyu Zhang

# Code #
* Link: https://github.com/lneef/PSEMolDyn_Group1
* Branch: main
* Revision:


# Run Instructions #
* for cluster: load cmake/3.14.5 and xerces/3.2.1, disable tests via -DTESTS=OFF
* Compiler: GNU 12.2.0, 11.3.0, 11.2.0, icpc 2021.4.0
* run mkdir build and cd build/
* run cmake ..
* run make or make all
* run ./MolSim --help for further information on usage or see section Usage
* run ctest for running all tests

# Doxygen #
* Per default cmake .. creates a target for generating a doxygen documentation
* you can disable this behaviour by passing DISABLE\_DOXYGEN=ON
* to create the doxygen documentation run make doc\_doxygen

# CMake #
* by setting SANITIZE=ON compilation with fsanitize=address is enables
* by setting BMARK=ON compilation in benchmark mode is enabled, i.e. all io in benchmark region is disabled
* by setting VECTOR=ON compilation with vectorization is enabled
* by setting DTESTS=OFF test target won't be build
* by setting FAST_I=ON and using the icpc compiler you can compile with high optimization using icpc
* by setting FAST_G=ON and using gcc you can compile with high optimization using gcc
* by setting OMP=ON compilation with openmp is activated for gcc
* by setting IOMP=ON compilation with openmp is activated for icpc

# Usage #
MolSim is a simulator for simulating the movement of planets and molecules.

The following command line options are supported:

    -e      Expects an argument of type double. Using this option the end time of the simulation can be set.
            The default end time is 1000

    -t      Expects an argument of type double. Using this option the time step (delta\_t) of the simulation can be set.
            The default time step is 0.014

    -l     Expects as argument either off, info, debug, trace or error. Using this option the global log_level can be set

    --planet Usage: --planet path-to-file. Using this command line option particles are read from the specified file and
             the planet simulation is started.

    --cub   Using this option the molecular simulation is started. Particle are generated from cuboid which are either read
            from a file or from the command line. In the former case --cub=path-to-file has to be used and in the latter
            --cub.

    --xml   Usage: --xml path-to-xml-file. Using this option all data for the simualtion can be read from an xml file.
            -e and -t are ignored by this option.

    --check   Usage: --check output-file-name. Using this option to add a checkpointing-file(.txt) witm the name   output-file-name

    --help  This page is printed.

    You must use the options --planet or --cub if you want to run a simulation. 

# Logging #
* One Logger for MolSim(MolSimLogger) which writes to logs.txt
* One Logger for tests(MolSimLogger_test) with loglevel debug which writes to logs_test.txt

# Input #
* for cuboids see input/eingabe-cuboid.txt
* for planets see input/eingabe-sonne.txt
* for cuboids via xml see input/input-cuboid.xml
* for spheres via xml see input/sphere.xml and the specification in sphere.txt
* for fluids via xml see input/instability_small.xml and input/instability_big.xml
* for trying checkpointing via xml see input/checkpoint.xml and then the input/from_checkpoint.xml with the checkpoint parameters in input/checkpoint.txt
* for our profiling setup via xml see input/cluster.xml
* two parallelization strategies (generic and tasking), use parallel_mode to specify which one you want to use(default generic)

# Boundary #
* the default boundary condition is outflow
* if you want to use reflection, write reflecting for the respective boundary into the xml file

# XML #
* see molsim.xsd for the schema definition
* simulation has to be the first element since it is needed for initialization
* to get an idea have a look at the provided files in input/
* for temperature the target temperature must appear after the initial temperature
* an element named path expects a valid path to an input file (.txt, see sphere.txt and eingabe-cuboid.txt for schema)

# Statistics #
* if you want to use the statistics, you must make sure that the folders /output, /output/rdf and /output/plot exist
* add the statistics-parameters to your xml-input file
* the statistics only work with the xml-input
* after you have run your simulation, you should run /src/plotStatistics.py
* you will find the plots as png files in /output/plot

# Plots #
the RDF- and Difussion graphs of both Cooling and Supercooling can be found at:
https://syncandshare.lrz.de/getlink/fiHf5MvzcJGzSjPM4iHRuf/

# Simulation Runs #
Videos of runs of the simulation be found at:
https://syncandshare.lrz.de/getlink/fiHf5MvzcJGzSjPM4iHRuf/

# Contest #
* the parallel_mode is tasking
* use icpc as described in Run Instructions and run cmake with -DBMARK=ON -DIOMP=ON -DTESTS=OFF (2D: 4531 ms and 2,2 mmups, 3D: 318112ms, 314355 mups)
* use gcc as described above and run cmake wih -DBMARK=ON -DOMP=ON (2D: 4603ms and 2,1 mmups, 3D 189764ms and 526970 mups)
