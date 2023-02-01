#include <iostream>
#include "MembraneForce.h"
#include "container/ParticleContainer.h"


MembraneForce::~MembraneForce() = default;

void MembraneForce::membraneForce(std::shared_ptr<Container> &particles, double f_up, double g, std::unique_ptr<Force> &force, double time){
    if (time <= 150) {
                particles->apply([g, f_up](Particle &p) {
                    if ((p.getIndex()[0] == 17 && p.getIndex()[1] == 24) ||
                        (p.getIndex()[0] == 17 && p.getIndex()[1] == 25) ||
                        (p.getIndex()[0] == 18 && p.getIndex()[1] == 24) ||
                        (p.getIndex()[0] == 18 && p.getIndex()[1] == 25)) {
                        p.updateF({0, 0, p.getM() * g + f_up});
                        // std::cout << p.getIndex().at(0) << " " << p.getIndex().at(1) << std::endl;
                    } else {
                        p.updateF({0, 0, p.getM() * g});
                    }
                });

            } else {
                particles->apply([g](Particle &p) {
                        p.updateF({0, 0,  p.getM() * g});

                });
            }

            force->calculateF(particles);
}


void MembraneForce::calculateF(std::shared_ptr<Container>& particles) {
    auto r0_arg = r0;
    auto k_arg = k;
    particles->applyF([r0_arg, k_arg](Particle& p1, Particle& p2) {
        calculateF(p1, p2, r0_arg, k_arg);
        });
}


MembraneForce::MembraneForce(double r0_arg, double k_arg) {
    r0 = r0_arg;
    k = k_arg;
}

void MembraneForce::calculateF(Particle& p1, Particle& p2, double r0, double k) {
    double sigma, epsilon;
    if (p1.getType() == p2.getType()) {
        sigma = p1.getSigma();
        epsilon = p1.getEpsilon();
    }
    else {
        sigma = (p1.getSigma() + p2.getSigma()) / 2.;
        epsilon = std::sqrt(p1.getEpsilon() * p2.getEpsilon());
    }
    std::array<double, 3> xij = p1.getX() - p2.getX();
    double sum = xij[0] * xij[0] + xij[1] * xij[1] + xij[2] * xij[2];

    if (p1.ifDirectNeighbor(p2)) {
        double norm = std::sqrt(sum);
        double scalar = k * (r0 / norm - 1);
        std::array<double, 3> newF = scalar * xij;
        p1.addToF(newF);
        p2.subtractFromF(newF);
    }
    else if (p1.ifDiagonalNeighbor(p2)) {
        double norm = std::sqrt(sum);
        double scalar = k * (1.414 * r0 / norm - 1);
        std::array<double, 3> newF = scalar * xij;
        p1.addToF(newF);
        p2.subtractFromF(newF);
    }
    else {
        if (sum < 1.26 * sigma * sigma) {
            double pow_6 = pow((sigma * sigma) / sum, 3);
            double scalar = ((-24 * epsilon) / sum) * (pow_6 - 2 * pow(pow_6, 2));
            std::array<double, 3> newF = scalar * xij;
            p1.addToF(newF);
            p2.subtractFromF(newF);
        }
    }
}
