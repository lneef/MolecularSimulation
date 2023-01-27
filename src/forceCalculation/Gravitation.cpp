

#include "Gravitation.h"

//calculation new force
void Gravitation::calculateF(std::shared_ptr<Container> &particles) {

    particles->apply([](Particle &p) {
        p.updateF({});
    });

    particles->applyF([](Particle &p1, Particle &p2) {
                          std::array<double, 3> xij = p1.getX() - p2.getX();
                          std::array<double, 3> xji = p2.getX() - p1.getX();
                          double norm = ArrayUtils::L2Norm(xij);
                          std::array<double, 3> newF{};
                          double scalar = p1.getM() * p2.getM() / pow(norm, 3);

                          newF = scalar * xji;
                          p1.addToF(newF);
                          p2.subtractFromF(newF);
                      }
    );
}

Gravitation::~Gravitation() = default;
