//
// Created by dominik on 11.01.23.
//

#include "Statistics.h"

Statistics::Statistics(int i_begin, int i_end, double delta_distance) {
    this->i_rdf_begin = i_begin;
    this->i_rdf_end = i_end;
    this->delta_r = delta_distance;
}

void Statistics::calcDiffusion(double delta_t) {
    double numerator = 0.;
    std::array<double, 3> dom_size = particles->getDomain();
    particles->apply([&numerator, &dom_size, &delta_t](Particle &p) {
        std::array<double, 3> x_diff{};
        std::array<double, 3> v = p.getV();
        std::array<double, 3> x = p.getX();
        std::array<double, 3> old_x = p.getOldX();
        std::array<double, 3> tempF = p.getF();
        double m = p.getM();
        for (int i = 0; i < 3; i++) {
            double checkValue = v[i] + delta_t * tempF[i] / (2 * m);
            if (x[i] < old_x[i] && checkValue > 0) {
                x_diff[i] = dom_size[i] - (old_x[i] - x[i]);
            } else if (x[i] > old_x[i] && checkValue < 0) {
                x_diff[i] = dom_size[i] - (x[i] - old_x[i]);
            } else {
                x_diff[i] = x[i] - old_x[i];
            }
        }
        double sum = x_diff[0] * x_diff[0] + x_diff[1] * x_diff[1] + x_diff[2] * x_diff[2];
        numerator = numerator + sum;
    });
    if (particles->size() > 0) {
        diffusion.push_back(numerator / particles->size());
    }
}

void Statistics::calcRDF() {
    int array_size = i_rdf_end - i_rdf_begin + 1;
    std::vector<double> quantities(array_size, 0.);
    std::vector<std::array<double, 3>> positions;
    particles->apply([&positions](Particle &p) {
        positions.push_back(p.getX());
    });
    std::vector<double> distances;
    for (size_t i = 0; i < positions.size(); ++i) {
        for (size_t j = 0; j < positions.size(); ++j) {
            if(i!=j){
                distances.push_back(sqrt(pow(positions[j][0] - positions[i][0], 2) + pow(positions[j][1] - positions[i][1], 2) + pow(positions[j][2] - positions[i][2], 2)));
            }
        }
    }
    for (
            int i = 0;
            i < array_size;
            i++) {
        for (
            double d
                : distances) {
            if ((d > (i + i_rdf_begin)) && (d <= (i + i_rdf_begin + delta_r))) {
                quantities[i] = quantities[i] + 1;
            }
        }
    }

    std::vector<double> loc_densities;
    for (
            int i = 0;
            i < array_size;
            i++) {
        loc_densities.
                push_back(
                (quantities[i]
                 / 2) /
                (((4 * M_PI) / (3)) * (
                        pow(i
                            + i_rdf_begin + delta_r, 3) -
                        pow(i
                            + i_rdf_begin, 3))));
    }
    rdf.
            push_back(loc_densities);
}

void Statistics::writeDiffusion() {
    std::ifstream checkFile;
    checkFile.open("../output/diffusion.csv");
    if (checkFile) {
        remove("../output/diffusion.csv");
    }
    std::ofstream diffFile;
    diffFile.open("../output/diffusion.csv");
    diffFile << "timestep, var(t)" << std::endl;
    for (int i = 0; i < diffusion.size(); i++) {
        diffFile << i * n_statistics << "," << diffusion[i] << std::endl;
    }
    diffFile.close();
}

void Statistics::writeRDF() {
    for (int i = 0; i < rdf.size(); i++) {
        std::ofstream rdfFile;
        std::string path = "../output/rdf/time_" + std::to_string(i * n_statistics) + ".csv";
        rdfFile.open(path);
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