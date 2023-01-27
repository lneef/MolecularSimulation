#include <gmock/gmock-matchers.h>
#include "forceCalculation/SLennardJones.h"
#include "gtest/gtest.h"
#include <iostream>

/**
 * @brief test whether calculateF of Smoothed Lennard Jones works a expected and old_f is updated as expected wenn dij<=rl
 */
TEST(SLennardJonesTest, SLJTest1) {
    std::shared_ptr<Container> par = std::make_shared<ParticleContainer>();
    par->addParticle(Particle({ 0.,0.,0. }, { 0.,0., 0. }, 1.));
    par->addParticle(Particle({ 1.0, 0. , 0. }, { 1.5, 0., 0. }, 1.));

    std::array<double, 3> oldf{ 1.,1., 1. };
    auto& p = dynamic_cast<ParticleContainer&>(*par);
    p[0].setF(oldf);

    SLennardJones slj{ 2.3,1.9 };

    slj.calculateF(par);
    const std::array<double, 3> f1{ -120., 0. ,0. };
    const std::array<double, 3> f2{ 120, 0., 0. };

    EXPECT_THAT(p[0].getF(), testing::Pointwise(testing::DoubleEq(), f1));
    EXPECT_THAT(p[1].getF(), testing::Pointwise(testing::DoubleEq(), f2));
    EXPECT_THAT(p[0].getOldF(), testing::Pointwise(testing::DoubleEq(), oldf));


}


/**
 * @brief test whether calculateF of Smoothed Lennard Jones works a expected and old_f is updated as expected wenn rl< dij<rc
 */
TEST(SLennardJonesTest, SLJTest2) {
    std::shared_ptr<Container> par = std::make_shared<ParticleContainer>();
    par->addParticle(Particle({ 0.,0.,0. }, { 0.,0., 0. }, 1.));
    par->addParticle(Particle({ 1.0, 0. , 0. }, { 1.5, 0., 0. }, 1.));

    std::array<double, 3> oldf{ 1.,1., 1. };
    auto& p = dynamic_cast<ParticleContainer&>(*par);
    p[0].setF(oldf);

    SLennardJones slj{ 1.5,0.5 };

    slj.calculateF(par);

    // p[0].setF({0.1,0,0});
    const std::array<double, 3> f1{ -60., 0. ,0. };
    const std::array<double, 3> f2{ 60, 0., 0. };

    EXPECT_THAT(p[0].getF(), testing::Pointwise(testing::DoubleEq(), f1));
    EXPECT_THAT(p[1].getF(), testing::Pointwise(testing::DoubleEq(), f2));
    EXPECT_THAT(p[0].getOldF(), testing::Pointwise(testing::DoubleEq(), oldf));


}

/**
 * @brief test whether calculateF of Smoothed Lennard Jones works a expected and old_f is updated as expected wenn dij>rc
 */
TEST(SLennardJonesTest, SLJTest3) {
    std::shared_ptr<Container> par = std::make_shared<ParticleContainer>();
    par->addParticle(Particle({ 0.,0.,0. }, { 0.,0., 0. }, 1.));
    par->addParticle(Particle({ 1.0, 0. , 0. }, { 1.5, 0., 0. }, 1.));

    std::array<double, 3> oldf{ 1.,1., 1. };
    auto& p = dynamic_cast<ParticleContainer&>(*par);
    p[0].setF(oldf);

    SLennardJones slj{ 0.9,0.5 };

    slj.calculateF(par);

    // p[0].setF({0.1,0,0});
    const std::array<double, 3> f1{ 0., 0. ,0. };
    const std::array<double, 3> f2{ 0, 0., 0. };

    EXPECT_THAT(p[0].getF(), testing::Pointwise(testing::DoubleEq(), f1));
    EXPECT_THAT(p[1].getF(), testing::Pointwise(testing::DoubleEq(), f2));
    EXPECT_THAT(p[0].getOldF(), testing::Pointwise(testing::DoubleEq(), oldf));


}