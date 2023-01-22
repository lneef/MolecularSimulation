//
// Created by dominik on 11.01.23.
//

#pragma once

#include <iostream>
#include <vector>
#include <array>
#include "container/LinkedCellContainer.h"
#include <cmath>

class Statistics {
private:
    std::vector<double> diffusion;
    std::vector<std::vector<double>> rdf;
    int i_rdf_begin;
    int i_rdf_end;
    double delta_r;
    std::shared_ptr<LinkedCellDataStructure> particles;
public:
    Statistics(int i_begin, int i_end, double delta_distance);

    void calcDiffusion();

    void calcRDF();

    void setParticles(std::shared_ptr<LinkedCellDataStructure> particles_arg);

    std::vector<std::vector<double>> getRdf() const;

    std::vector<double> getDiffusion() const;

};