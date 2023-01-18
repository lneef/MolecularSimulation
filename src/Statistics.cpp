//
// Created by dominik on 11.01.23.
//

#include "Statistics.h"

Statistics::Statistics(int i_begin, int i_end, double delta_distance) {
    this->i_rdf_begin = i_begin;
    this->i_rdf_end = i_end;
    this->delta_r = delta_distance;
}

void Statistics::calcDiffusion() {
    double numerator = 0.;
    particles->apply([&numerator](Particle &p) {
        //ToDo: perodic boundary
        std::array<double, 3> x_diff = p.getX() - p.getOldX();
        double sum = x_diff[0] * x_diff[0] + x_diff[1] * x_diff[1] + x_diff[2] * x_diff[2];
        numerator = numerator + sum;
    });
    diffusion.push_back(numerator / particles->size());
}

void Statistics::calcRDF() {
    int array_size = i_rdf_end - i_rdf_begin + 1;
    double quantities[array_size];
    std::vector<std::array<double, 3>> positions;
    //ToDo Apply pairs of particles
    particles->apply([&positions](Particle &p) {
        positions.push_back(p.getX());
    });
    std::vector<double> distances;
    for (std::array<double, 3> x1: positions) {
        for (std::array<double, 3> x2: positions) {
            if (x1 != x2) {
                distances.push_back(sqrt(pow(x2[0] - x1[0], 2) + pow(x2[1] - x1[1], 2) + pow(x2[2] - x1[2], 2)));
            }
        }
    }
    for (int i = 0; i < array_size; i++) {
        for (double d: distances) {
            if (d > i + i_rdf_begin && d <= i + i_rdf_begin + delta_r) {
                quantities[i]++;
            }
        }
    }
    std::vector<double> loc_densities;
    for(int i = 0; i < array_size; i++){
        loc_densities.push_back((quantities[i])/(((4*M_PI)/(3))*(pow(i+i_rdf_begin+delta_r,3)-pow(i+i_rdf_begin,3))));
    }
    rdf.push_back(loc_densities);
}

void Statistics::setParticles(std::shared_ptr<Container> particles_arg) {
    particles = particles_arg;
}
