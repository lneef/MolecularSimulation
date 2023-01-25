#include "gtest/gtest.h"
#include "Statistics.h"
#include "container/LinkedCell3D.h"

TEST(RDFTest, StatisticsTest) {
    std::shared_ptr<LinkedCell3D> lc = std::make_shared<LinkedCell3D>();
    std::array<double, 3> domain = {10., 10., 1.};
    lc->setDomain(domain);

    std::array<double, 3> x;
    std::array<double, 3> v;
    std::array<double, 3> old_x;
    double m = 1.;

    /*lc->addParticle(Particle({5.,5.,0.},{0.,0.,0.},1.));
    lc->addParticle(Particle({7.,5.,0.},{0.,0.,0.},1.));
    lc->addParticle(Particle({3.,5.,0.},{0.,0.,0.},1.));

    Statistics s(1,3,1.);
    s.setParticles(lc);
    s.calcRDF();
    std::vector<std::vector<double>> rdf = s.getRdf();
    EXPECT_EQ(round(rdf[0][0] * 1000), 68.);
    EXPECT_EQ(round(rdf[0][2] * 1000), 6.);*/
}