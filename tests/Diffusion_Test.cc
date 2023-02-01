#include "gtest/gtest.h"
#include "Statistics.h"
#include "container/LinkedCell3D.h"

TEST(DiffusionTest, StatisticsTest) {
    std::shared_ptr<LinkedCell3D> lc = std::make_shared<LinkedCell3D>();
    std::array<double,3> domain = {10.,10., 1.};
    lc->setSize(1,domain);

    std::array<double, 3> x;
    std::array<double, 3> v;
    std::array<double, 3> old_x;
    double m = 1.;

    //periodic right
    x = {2., 5., 0.};
    old_x = {8., 5., 0.};
    v = {3., 0., 0.};
    Particle p1(x, v, m);
    p1.setOldX(old_x);
    lc->addParticle(p1);

    //periodic left
    x = {8., 5., 0.};
    old_x = {2., 5., 0.};
    v = {-3., 0., 0.};
    Particle p2(x, v, m);
    p2.setOldX(old_x);
    lc->addParticle(p2);

    //periodic corner
    x = {2., 2., 0.};
    old_x = {8., 8., 0.};
    v = {3., 3., 0.};
    Particle p3(x, v, m);
    p3.setOldX(old_x);
    lc->addParticle(p3);

    //normal x
    x = {8., 5., 0.};
    old_x = {2., 5., 0.};
    v = {3., 0., 0.};
    Particle p4(x, v, m);
    p4.setOldX(old_x);
    lc->addParticle(p4);

    Statistics s(1,3,1.);
    s.setParticles(lc);
    s.calcDiffusion();
    std::vector<double> diffusion = s.getDiffusion();
    EXPECT_EQ(diffusion[0], 45.);
}