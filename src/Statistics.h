//
// Created by dominik on 11.01.23.
//

#pragma once

#include <iostream>
#include <vector>
#include <array>
#include "container/LinkedCellContainer.h"

class Statistics {
private:
    std::vector<double> diffusion;
    std::vector<std::vector<int>> rdf;
    int i_rdf_begin;
    int i_rdf_end;
    double delta_r;
    std::shared_ptr<Container> particles;
public:
    Statistics(int i_begin, int i_end, double delta_distance, std::shared_ptr<Container> &particles);
    void calcDiffusion();
    void calcRDF();
};