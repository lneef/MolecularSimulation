//
// Created by lukas on 08.11.22.
//
#pragma once

#include <memory>
#include <iostream>


#include "outputWriter/FileWriter.h"
#include "forceCalculation/Force.h"
#include "container/LinkedCellContainer.h"
#include "tempCalculation/Thermostat.h"
#include "tempCalculation/Temperature.h"
#include <chrono>

/**
 * @brief The class Simulation provides the functionality for running the simulation. It can be tailored to a specific
 * use case by initializing the fields accordingly
 */
class Simulation {
    /**
     * @brief Container containing the particles for the simulation
     */
    std::shared_ptr<Container> particles;
    /**
     * @brief start time of the simulation
     */
    const double start_time = 0.;

    /**
     * @brief time step after which values are updated
     */
    double delta_t;

    /**
     * @brief time at which the simulation ends
     */
    double end_time;

    /**
     * @brief FileWriter to print the VTK files
     */
    std::unique_ptr<outputWriter::FileWriter> writer;

    /**
     * @brief Force to calculate the force between particles
     */
    std::unique_ptr<Force> force;

    /**
     * @brief filename of the output files (default: MD_vtk)
     */
    std::string out_name = "MD_vtk";

    /**
     * @brief output frequency (default: every 10 iterations)
     */
    int out_frequency = 10;

    std::unique_ptr<Temperature> temperature;
    std::shared_ptr<Thermostat> thermostat;
    int n_thermostat;

public:
    /**
     * @brief Calculates next position for each particle in particles
     */
    void calculateX();

    /**
     * @brief Calculates next velocity for each particle in particles
     */
    void calculateV();

    /**
     * @brief Function to start and run the simulation
     */
    void run();

    /**
     * @brief Constructor to initialize the simulation based on the use case.
     * @param particles Container to store the involved particles
     * @param delta_t time step for which the position, velocity and force is recalculated
     * @param end_time time at which the simulation ends
     * @param writer object to specify the output behaviour of the simulation
     * @param force object to specify the type of force calculation used during the simulation
     */

    Simulation(std::shared_ptr<Container> &particles, double delta_t, double end_time,
               std::unique_ptr<outputWriter::FileWriter> &writer, std::unique_ptr<Force> &force);

    /**
     * @brief constructor of Simulation which only initializes delta_t and end_time
     * @param delta_t_arg duration of a time step
     * @param end_time_arg end time of the simulation
     */
    explicit Simulation(double delta_t_arg = 2, double end_time_arg = 0.0002);

    /**
     * @brief setter for delta_t
     * @param delta_t_arg duration of a timestep of the simulation
     */
    void setDeltaT(double delta_t_arg);

    /**
     * @brief setter for end_time
     * @param end_time_arg end time of the simulation
     */
    void setEndTime(double end_time_arg);

    /**
     * @brief setter for force
     * @param force_arg method for calculating the force effective on particles
     */
    void setForce(std::unique_ptr<Force> &force_arg);

    /**
     * @brief setter for particles
     * @param particles_arg ParticleContainer containing the particles
     */
    void setParticle(std::shared_ptr<ParticleContainer> &particles_arg);

    /**
     * @brief setter for particles
     * @param particles_arg LinkedCellContainer containing the particles
     */
    void setParticle(std::shared_ptr<LinkedCellContainer> &particles_arg);

    /**
     * @brief setter for out_name
     * @param out_name_arg base name for output file
     */
    void setOut_name(const std::string &out_name_arg);

    /**
     * @brief setter for output frequency
     * @param out_frequency_arg output frequency of the simulation
     */
    void setOut_frequency(int out_frequency_arg);

    /**
     * @brief setter for output writer
     * @param writer_arg file writer for writing output to a file
     */
    void setWriter(std::unique_ptr<outputWriter::FileWriter> &writer_arg);

    void setN_thermostat(int n_thermostat);
    void setThermostat(std::shared_ptr<Thermostat> &thermostat);
    void setTemperature(std::unique_ptr<Temperature> &temperature);
    const std::shared_ptr<Thermostat> &getThermostat() const;
};



