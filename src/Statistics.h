//
// Created by dominik on 11.01.23.
//

#pragma once

#include <iostream>
#include <vector>
#include <array>
#include "container/LinkedCellContainer.h"
#include <cmath>
#include <fstream>
#include <filesystem>
#include "utils/ArrayUtils.h"

class Statistics {
private:
    /**
     * @brief Vector of the calculated diffusion values
     */
    std::vector<double> diffusion;
    /**
     * @brief Vector of vectors for each timestep with the calculated rdf values
     */
    std::vector<std::vector<double>> rdf;
    /**
     * @brief First value of the rdf distance intervals
     */
    int i_rdf_begin;
    /**
     * @brief Last value of the rdf distance intervals
     */
    int i_rdf_end;
    /**
     * @brief Delta r, which defines the length of the rdf distance intervals
     */
    double delta_r;
    /**
     * @brief Shared pointer of the particles in the simulation
     */
    std::shared_ptr<LinkedCellDataStructure> particles;
    /**
     * @brief Time step after which statistics will be calculated
     */
    int n_statistics;
public:
    /**
     * @brief Constructor of the class Statistics
     * @param First value of the rdf distance intervals
     * @param Last value of the rdf distance intervals
     * @param Delta r, which defines the length of the rdf distance intervals
     */
    Statistics(int i_begin, int i_end, double delta_distance);

    /**
     * @brief Calculates the diffusion values and stores them in the diffusion vector
     */
    void calcDiffusion();
    /**
     * @brief Calculates the rdf values and stores them in a vector for this calculation iteration, which is then stored in the rdf vector
     */
    void calcRDF();
    /**
     * @brief Setter for particles
     * @param Shared pointer of the particles in the simulation
     */
    void setParticles(std::shared_ptr<LinkedCellDataStructure> particles_arg);
    /**
     * @brief Getter for the rdf vector
     * @return Rdf vector
     */
    std::vector<std::vector<double>> getRdf() const;
    /**
     * @brief Getter for the diffusion vector
     * @return Diffusion vector
     */
    std::vector<double> getDiffusion() const;
    /**
     * @brief Writes the values of the diffusion vector in a csv file
     */
    void writeDiffusion();
    /**
     * @brief Writes the values of the rdf vector in csv files
     */
    void writeRDF();

    /**
     * @brief Setter for n_statistics
     * @param Time step after which statistics will be calculated
     */
    void setN_statistics(int n_arg);
};