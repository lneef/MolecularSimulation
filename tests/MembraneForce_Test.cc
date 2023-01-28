#include <gmock/gmock-matchers.h>
#include "forceCalculation/MembraneForce.h"
#include "gtest/gtest.h"
#include <iostream>

/**
 * @brief test whether calculateF of MembraneForce works a expected and old_f is updated as expected wenn they are direct neighbors
 */
TEST(MembraneForceTest, MFTest1) {
    std::shared_ptr<Container> par = std::make_shared<ParticleContainer>();
    par->addParticle(Particle({ 0.,0.,0. }, { 0.,0., 0. }, 1.));
    par->addParticle(Particle({ 1.0, 0. , 0. }, { 1.5, 0., 0. }, 1.));

    std::array<double, 3> oldf{ 1.,1., 1. };
    auto& p = dynamic_cast<ParticleContainer&>(*par);
    p[0].setF(oldf);
    p[0].setIndex({ 0,0 });
    p[1].setIndex({ 0,1 });
    MembraneForce mf{ 2.3,300 };

    mf.calculateF(par);
    const std::array<double, 3> f1{ -390., 0. ,0. };
    const std::array<double, 3> f2{ 390, 0., 0. };

    EXPECT_THAT(p[0].getF(), testing::Pointwise(testing::DoubleEq(), f1));
    EXPECT_THAT(p[1].getF(), testing::Pointwise(testing::DoubleEq(), f2));
    EXPECT_THAT(p[0].getOldF(), testing::Pointwise(testing::DoubleEq(), oldf));


}

/**
 * @brief test whether calculateF of MembraneForce works a expected and old_f is updated as expected wenn they are diagonal neighbors
 */
TEST(MembraneForceTest, MFTest2) {
    std::shared_ptr<Container> par = std::make_shared<ParticleContainer>();
    par->addParticle(Particle({ 0.,0.,0. }, { 0.,0., 0. }, 1.));
    par->addParticle(Particle({ 1.0, 0. , 0. }, { 1.5, 0., 0. }, 1.));

    std::array<double, 3> oldf{ 1.,1., 1. };
    auto& p = dynamic_cast<ParticleContainer&>(*par);
    p[0].setF(oldf);
    p[0].setIndex({ 0,0 });
    p[1].setIndex({ 1,1 });
    MembraneForce mf{ 2.3,300 };

    mf.calculateF(par);
    const std::array<double, 3> f1{ -675.66, 0. ,0. };
    const std::array<double, 3> f2{ 675.66, 0., 0. };

    EXPECT_THAT(p[0].getF(), testing::Pointwise(testing::DoubleEq(), f1));
    EXPECT_THAT(p[1].getF(), testing::Pointwise(testing::DoubleEq(), f2));
    EXPECT_THAT(p[0].getOldF(), testing::Pointwise(testing::DoubleEq(), oldf));
}


/**
 * @brief test whether calculateF of MembraneForce works a expected and old_f is updated as expected wenn they are not neighbors, but sill has forces
 */
TEST(MembraneForceTest, MFTest3) {
    std::shared_ptr<Container> par = std::make_shared<ParticleContainer>();
    par->addParticle(Particle({ 0.,0.,0. }, { 0.,0., 0. }, 1.));
    par->addParticle(Particle({ 1.0, 0. , 0. }, { 1.5, 0., 0. }, 1.));

    std::array<double, 3> oldf{ 1.,1., 1. };
    auto& p = dynamic_cast<ParticleContainer&>(*par);
    p[0].setF(oldf);
    p[0].setIndex({ 0,0 });
    p[1].setIndex({ 1,2 });
    MembraneForce mf{ 2.3,300 };

    mf.calculateF(par);
    const std::array<double, 3> f1{ -120, 0. ,0. };
    const std::array<double, 3> f2{ 120, 0., 0. };

    EXPECT_THAT(p[0].getF(), testing::Pointwise(testing::DoubleEq(), f1));
    EXPECT_THAT(p[1].getF(), testing::Pointwise(testing::DoubleEq(), f2));
    EXPECT_THAT(p[0].getOldF(), testing::Pointwise(testing::DoubleEq(), oldf));
}

/**
 * @brief test whether calculateF of MembraneForce works a expected and old_f is updated as expected wenn they are not neighbors, but has no forces
 */
TEST(MembraneForceTest, MFTest4) {
    std::shared_ptr<Container> par = std::make_shared<ParticleContainer>();
    par->addParticle(Particle({ 0.,0.,0. }, { 0.,0., 0. }, 1.));
    par->addParticle(Particle({ 1.27, 0. , 0. }, { 1.5, 0., 0. }, 1.));

    std::array<double, 3> oldf{ 1.,1., 1. };
    auto& p = dynamic_cast<ParticleContainer&>(*par);
    p[0].setF(oldf);
    p[0].setIndex({ 0,0 });
    p[1].setIndex({ 1,2 });
    MembraneForce mf{ 2.3,300 };

    mf.calculateF(par);
    const std::array<double, 3> f1{ 0, 0. ,0. };
    const std::array<double, 3> f2{ 0, 0., 0. };

    EXPECT_THAT(p[0].getF(), testing::Pointwise(testing::DoubleEq(), f1));
    EXPECT_THAT(p[1].getF(), testing::Pointwise(testing::DoubleEq(), f2));
    EXPECT_THAT(p[0].getOldF(), testing::Pointwise(testing::DoubleEq(), oldf));
}
