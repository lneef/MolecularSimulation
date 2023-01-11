#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include "inputReader/xmlReader/XmlReader.h"

TEST(MembraneGeneratorTest, GenerationTest) {
    ParticleGenerator<ParticleContainer> cg;
    std::shared_ptr<ParticleContainer> par = std::make_shared<ParticleContainer>();
    std::array<double, 3> x = {0., 0., 0.};
    std::array<int, 3> n = {25, 25, 1};
    double h = 1.1;
    double m = 1.;
    std::array<double, 3> v = {0., 0., 0.};

    cg.generateMembraneBrownian(par, x, n, v, h, m, 1, 1, 5, 3, 0.8);
    std::array<double, 3> pos{0, 1.1, 0};
    std::array<double, 3> zero_force{};
    std::array<double, 3> fz_force{0, 0, 0.8};
    std::array<int, 2> index{0, 1};
    EXPECT_THAT((*par)[1].getX(), testing::Pointwise(testing::DoubleEq(), pos));
    EXPECT_THAT((*par)[1].getF(), testing::Pointwise(testing::DoubleEq(), zero_force));
    EXPECT_EQ((*par)[1].getIndex(), index);

    pos = {18.7, 26.4, 0};
    index = {17, 24};
    EXPECT_THAT((*par)[449].getX(), testing::Pointwise(testing::DoubleEq(), pos));
    EXPECT_THAT((*par)[449].getF(), testing::Pointwise(testing::DoubleEq(), fz_force));
    EXPECT_EQ((*par)[449].getIndex(), index);

}