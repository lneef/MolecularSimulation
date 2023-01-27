#include <iostream>
#include "SLennardJones.h"
#include "container/ParticleContainer.h"


SLennardJones::~SLennardJones() = default;

//calculation new force with harmonic potential

void SLennardJones::calculateF(std::shared_ptr<Container>& particles) {

    particles->apply([](Particle& p) {
        p.updateF({});
        });

    auto rc_arg = rc;
    auto rl_arg = rl;
    particles->applyF([rc_arg, rl_arg](Particle& p1, Particle& p2) {
        calculateF(p1, p2, rc_arg, rl_arg);
        });
}


SLennardJones::SLennardJones(double rc_arg, double rl_arg) {
    rc = rc_arg;
    rl = rl_arg;
}

void SLennardJones::calculateF(Particle& p1, Particle& p2, double rc, double rl) {

    std::array<double, 3> xij = p1.getX() - p2.getX();
    double sum = xij[0] * xij[0] + xij[1] * xij[1] + xij[2] * xij[2];
    double norm = std::sqrt(sum);

    if (norm >= rc) {
        return;
    }

    double sigma, epsilon;
    if (p1.getType() == p2.getType()) {
        sigma = p1.getSigma();
        epsilon = p1.getEpsilon();
    }
    else {
        sigma = (p1.getSigma() + p2.getSigma()) / 2.;
        epsilon = std::sqrt(p1.getEpsilon() * p2.getEpsilon());
    }

    if (norm <= rl) {
        double pow_6 = pow((sigma * sigma) / sum, 3);
        double scalar = ((-24 * epsilon) / sum) * (pow_6 - 2 * pow(pow_6, 2));
        std::array<double, 3> newF = scalar * xij;
        p1.addToF(newF);
        p2.subtractFromF(newF);
    }
    else {
        double sigma_6 = pow(sigma, 6);
        double dij_6 = pow(sum, 3);

        double scalar1 = 24 * sigma_6 * epsilon / ((dij_6 * dij_6 * sum) * pow((rc - rl), 3)) * (rc - norm);
        double scalar2 = (rc + norm - 3 * rl) * (rc * (2 * sigma_6 - dij_6)) + norm * (5 * rl * sigma_6 - 2 * rl * dij_6 - 3 * sigma_6 * norm) + pow(sum, 4);
        std::array<double, 3> newF = scalar1 * scalar2 * xij;
        p1.addToF(newF);
        p2.subtractFromF(newF);
    }
}
