#include <iostream>
#include "LennardJones.h"
#include "container/ParticleContainer.h"
#include  <immintrin.h>

LennardJones::~LennardJones() = default;

//calculation new force with Lennard-Jones potential

void LennardJones::calculateF(std::shared_ptr<Container> &particles) {

    particles->apply([](Particle &p) {
        p.updateF({});
    });
    particles->applyF([](Particle &p1, Particle &p2) {
        calculateF(p1, p2);
    });

}


LennardJones::LennardJones() = default;

void LennardJones::calculateF(Particle &p1, Particle &p2) {

    double sigma, epsilon;

    if (p1.getType() == p2.getType()) {
        sigma = p1.getSigma();
        epsilon = p1.getEpsilon();
    } else {
        sigma = (p1.getSigma() + p2.getSigma()) / 2.;
        epsilon = std::sqrt(p1.getEpsilon() * p2.getEpsilon());
    }


#ifndef AVX
    std::array<double, 3> xij = p1.getX() - p2.getX();
    double sum = xij[0] * xij[0] + xij[1] * xij[1] + xij[2] * xij[2];
    double pow_6 = pow((sigma * sigma)/sum, 3);
    double scalar = ((-24 * epsilon) / sum) * (pow_6 - 2 * pow(pow_6, 2));
    std::array<double, 3> newF = scalar * xij;
    p1.addToF(newF);
    p2.subtractFromF(newF);
#endif

#ifdef AVX
    auto & pos1 = p1.getX();
    auto & pos2 = p2.getX();


    __m256d acc = _mm256_set_pd(pos1[2], pos1[1], pos1[0], 0);
    __m256d pp2 = _mm256_set_pd(pos2[2], pos2[1], pos2[0], 0);

    acc = _mm256_sub_pd(acc, pp2);


    __m256d norm_avx = _mm256_mul_pd(acc, acc);


    double norm = norm_avx[3] + norm_avx[2] + norm_avx[1];

    double pow_6 = pow(((sigma * sigma)/ norm), 3);
    double scalar = ((-24 * epsilon) / norm) * (pow_6 - 2 * pow(pow_6, 2));
    __m256d scalar_avx = _mm256_set_pd(scalar, scalar, scalar, scalar);
    acc = _mm256_mul_pd(acc, scalar_avx);
    std::array<double, 3> newF{acc[1], acc[2], acc[3]};
    p1.addToF(newF);
    p2.subtractFromF(newF);
#endif

}


