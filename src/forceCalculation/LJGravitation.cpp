

#include "LJGravitation.h"

LJGravitation::LJGravitation(double g_arg) {
    g = g_arg;
}

void LJGravitation::calculateF(std::shared_ptr<Container> &particles) {
    particles->applyPar([this](Particle& p){
        p.updateF({0, p.getM() * g, 0});
    });

    particles->applyF([](Particle& p1, Particle& p2){
        LennardJones::calculateF(p1, p2);
    });
}

LJGravitation::~LJGravitation() = default;