//
// Created by dominik on 11.01.23.
//

#include "Statistics.h"

Statistics::Statistics(int i_begin, int i_end, double delta_distance,std::shared_ptr<Container> &particles) {
    this->i_rdf_begin = i_begin;
    this->i_rdf_end = i_end;
    this->delta_r = delta_distance;
    this->particles = std::move(particles);
}

void Statistics::calcDiffusion(){
    double numerator = 0.;
    particles->apply([&numerator](Particle &p){
        //ToDo: perodic boundary
        std::array<double, 3> x_diff = p.getX() - p.getOldX();
        double sum = x_diff[0] * x_diff[0] + x_diff[1] * x_diff[1] + x_diff[2] * x_diff[2];
        numerator = numerator + sum;
    });
    diffusion.push_back(numerator/particles->size());
}

void Statistics::calcRDF() {
    /*int array_size = i_rdf_end - i_rdf_begin + 1;
    int quanities[array_size];


    std::vector<int> array_vec(quanities, quanities + sizeof quanities / sizeof quanities[0]);
    rdf.push_back(array_vec);*/
}
