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
    std::array<double, 3> dom_size = particles->getDomain();
    particles->apply([&numerator, &dom_size](Particle &p) {
        std::array<double, 3> x_diff;
        std::array<double, 3> v = p.getV();
        std::array<double, 3> x = p.getX();
        std::array<double, 3> old_x = p.getOldX();
        for (int i = 0; i < 3; i++) {
            if (x[i] < old_x[i] && v[i] > 0) {
                x_diff[i] = dom_size[i] - (old_x[i] - x[i]);
            } else if (x[i] > old_x[i] && v[i] < 0) {
                x_diff[i] = dom_size[i] - (x[i] - old_x[i]);
            } else {
                x_diff[i] = x[i] - old_x[i];
            }
        }
        double sum = x_diff[0] * x_diff[0] + x_diff[1] * x_diff[1] + x_diff[2] * x_diff[2];
        numerator = numerator + sum;
    });
    //if (particles->size() > 0) {
    diffusion.push_back(numerator / particles->size());
    //}
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
            if ((d > (i + i_rdf_begin)) && (d <= (i + i_rdf_begin + delta_r))) {
                quantities[i]++;
            }
        }
    }

    std::vector<double> loc_densities;
    for (int i = 0; i < array_size; i++) {
        loc_densities.push_back(
                (quantities[i] / 2) /
                (((4 * M_PI) / (3)) * (pow(i + i_rdf_begin + delta_r, 3) - pow(i + i_rdf_begin, 3))));
    }
    rdf.push_back(loc_densities);
}

void Statistics::writeDiffusion() {
    std::ifstream checkFile;
    checkFile.open("../output/diffusion");
    if (checkFile) {
        remove("../output/diffusion");
    }
    std::ofstream diffFile;
    diffFile.open("../output/diffusion");
    diffFile << "timestep, var(t)" << std::endl;
    for (int i = 0; i < diffusion.size(); i++) {
        diffFile << i * n_statistics << "," << diffusion[i] << std::endl;
    }
    diffFile.close();
}

void Statistics::writeRDF() {
    if (!std::filesystem::is_empty("../output/rdf")) {
        std::filesystem::remove_all("../output/rdf");
    }
    for (int i = 0; i < rdf.size(); i++) {
        std::ofstream rdfFile;
        rdfFile.open("../output/rdf/time_" + i * n_statistics);
        rdfFile << "distance, densities" << std::endl;
        for (int z = 0; z < rdf[i].size(); z++) {
            rdfFile << z + i_rdf_begin << "," << rdf[i][z] << std::endl;
        }
        rdfFile.close();
    }
}

void Statistics::setParticles(std::shared_ptr<LinkedCellDataStructure> particles_arg) {
    particles = particles_arg;
}

std::vector<std::vector<double>> Statistics::getRdf() const {
    return rdf;
}

std::vector<double> Statistics::getDiffusion() const {
    return diffusion;
}

void Statistics::setN_statistics(int n_arg) {
    n_statistics = n_arg;
}
