#include "gtest/gtest.h"
#include "container/LinkedCell3D.h"
#include "inputReader/ParticleGenerator.h"
#include "Simulation.h"


class LinkedCell3DTest : public ::testing::Test{
protected:
    std::shared_ptr<LinkedCell3D> lc;
    void SetUp() override{
        ParticleGenerator<LinkedCell3D> gen{};
        lc= std::make_shared<LinkedCell3D>();
        std::array<double, 3> dom{3, 3, 3};
        lc->setSize(1, dom);
        gen.generateCuboid(lc, {0.5, 0.5, 0.5},{3, 3, 3}, 1, 1, {0, 0, 1});
    }
    void TearDown() override {
        lc->clearBoundary();
    }
};
TEST_F(LinkedCell3DTest, AddTest){
    EXPECT_EQ(lc->size(), 27);
    EXPECT_EQ((*lc)[2][6].size(), 1);
    EXPECT_EQ((*lc)[2][12].size(), 1);
}

TEST_F(LinkedCell3DTest, BoundaryTest) {
    LinkedCell3D::addPeriodic(Boundary::BACK);
    LinkedCell3D::addPeriodic(Boundary::FRONT);
    lc->applyF([](Particle &p1, Particle &p2) {
        std::array<double, 3> add = {1., 0., 0.};
        p1.setF(p1.getF() + add);
        p2.setF(p2.getF() + add);
    });
    auto f = (*lc)[1][6].begin();
    auto & force = f->getF();
    auto f1 = (*lc)[3][6].begin();
    auto f2 = (*lc)[2][12].begin();
    auto & halo = (*lc)[0].getHalo();
    int sz = std::accumulate(halo.cbegin(), halo.cend(), 0, []( int cur, std::reference_wrapper<ParticleContainer> p){
       return cur + p.get().size();
    });
    EXPECT_EQ(sz, 16);
    EXPECT_EQ(force[0], 16);
    EXPECT_EQ(f1->getF()[0], 16);
    EXPECT_EQ(f2->getF()[0], 26);

}

TEST_F(LinkedCell3DTest, BoundaryTest1) {
    LinkedCell3D::addPeriodic(Boundary::BACK);
    LinkedCell3D::addPeriodic(Boundary::FRONT);
    LinkedCell3D::addPeriodic(Boundary::LEFT);
    LinkedCell3D::addPeriodic(Boundary::RIGHT);
    lc->applyF([](Particle &p1, Particle &p2) {
        std::array<double, 3> add = {1., 0., 0.};
        p1.setF(p1.getF() + add);
        p2.setF(p2.getF() + add);
    });
    auto f = (*lc)[1][6].begin();
    auto & force = f->getF();
    auto & halo = (*lc)[0].getHalo();
    EXPECT_EQ(force[0], 23);

}

TEST_F(LinkedCell3DTest, MoveTest){
    LinkedCell3D::addPeriodic(Boundary::BACK);
    LinkedCell3D::addPeriodic(Boundary::FRONT);
    Simulation sim(1, 1);
    std::shared_ptr<LinkedCellDataStructure> test = lc;
    sim.setParticle(test);
    sim.setDeltaT(1);
    sim.calculateX();

    EXPECT_EQ(lc->size(), 27);
    EXPECT_EQ((*lc)[1][11].size(), 1);

}