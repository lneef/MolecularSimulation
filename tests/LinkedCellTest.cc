
#include <gmock/gmock-matchers.h>
#include "gtest/gtest.h"
#include "container/LinkedCellContainer.h"
#include "inputReader/ParticleGenerator.h"
#include "Simulation.h"
#include "MolSimLogger.h"

/**
 * @brief Test class with uniform setUP and Teardown function
 */
class LinkedCellTest : public testing::Test {
protected:
    std::shared_ptr<LinkedCellContainer> test;

    void SetUp() override {
        ParticleGenerator<LinkedCellContainer> cub{};
        test = std::make_shared<LinkedCellContainer>();
        std::array<double, 3> domain{3., 3., 0.5};
        test->setDomain(domain);
        test->setSize(1., domain);
        cub.generateCuboid(test, {.55, .55, 0}, {3, 3, 1}, 1., 1.0, {0., 1, 0.}, 1);
    }

    void TearDown() override {
        test->clearBoundary();
    }
};

/**
 * @brief tests if particles are added to the right cells
 */
TEST_F(LinkedCellTest, AddTest) {

    std::vector<ParticleContainer> celllist = test->getCells();
    EXPECT_EQ(test->size(), 9);
    EXPECT_EQ(celllist[6].size(), 1);
    EXPECT_EQ(celllist[7].size(), 1);
    EXPECT_EQ(celllist[11].size(), 1);
    EXPECT_EQ(celllist[12].size(), 1);


}

/**
 * @brief test if the correct number of references are stored
 */
TEST_F(LinkedCellTest, BoundaryHalo_Test){
    EXPECT_EQ(test->getBoundary().size(), 12);
    EXPECT_EQ(test->getHalo().size(), 16);
}

/**
 * @brief tests if force calculation function is applied to all particles correctly
 */
TEST_F(LinkedCellTest, AppTest) {
    test->applyF([](Particle &p1, Particle &p2) {
        std::array<double, 3> add = {1., 0., 0.};
        p1.addToF(add);
        p2.addToF(add);
    });
    auto celllist = test->getCells();
    auto it1 = celllist[6].begin();
    auto it2 = celllist[11].begin();
    auto it4 = celllist[12].begin();
    auto it3 = celllist[17].begin();
    
    EXPECT_DOUBLE_EQ((*it1).getF()[0], 3);
    MolSimLogger::logInfo("{}", (*it1).getF()[0]);
    EXPECT_DOUBLE_EQ(it2->getF()[0], 5);
    EXPECT_DOUBLE_EQ((*it4).getF()[0], 8);
    EXPECT_DOUBLE_EQ(it3->getF()[0], 5);

}

/**
 * @brief tests if reflecting boundary condition is applied to all particles
 */
TEST_F(LinkedCellTest, ReflectingBoundary){
    std::array<double ,3 > arr{1., 0, 0};
    test->addReflecting(Boundary::LEFT,Reflecting(arr, 0));
    test->addReflecting(Boundary::RIGHT, Reflecting(arr, 3.));
    test->applyF([](Particle &p1, Particle &p2) {
        std::array<double, 3> add = {1., 0., 0.};
        p1.addToF(add);
        p2.addToF(add);
    });
    std::vector<ParticleContainer> celllist = test->getCells();
    auto it1 = celllist[6].begin();
    auto it2 = celllist[8].begin();
    auto it4 = celllist[11].begin();


    EXPECT_DOUBLE_EQ((*it1).getF()[0], 4);
    EXPECT_DOUBLE_EQ((*it4).getF()[0], 6);
    EXPECT_DOUBLE_EQ((*it2).getF()[0], 4);
}

/**
 * @brief test if particles leaving the domain are deleted
 */
TEST_F(LinkedCellTest, Outflow){
    Simulation sim(1, 1);
    sim.setParticle(test);
    sim.calculateX();
    auto celllist = test-> getCells();
    EXPECT_EQ(test->size(), 6);
    EXPECT_EQ(celllist[11].size(), 1);
    EXPECT_EQ(celllist[7].size(), 0);
    EXPECT_EQ(celllist[18].size(), 1);


}

TEST_F(LinkedCellTest, Perioidc){
    test->addPeriodic(Boundary::TOP);
    test->addPeriodic(Boundary::BOTTOM);
    Simulation sim(1, 1);

    sim.setParticle(test);
    sim.calculateX();
    auto celllist = test-> getCells();
    EXPECT_EQ(test->size(), 9);
    EXPECT_EQ(celllist[6].size(), 1);
    EXPECT_EQ(celllist[7].size(), 1);
    EXPECT_EQ(celllist[8].size(), 1);
    EXPECT_EQ(celllist[21].size(), 0);


}

TEST_F(LinkedCellTest, PeriodicTest){
    test->addPeriodic(Boundary::LEFT);
    test->addPeriodic(Boundary::RIGHT);
    test->applyF([](Particle &p1, Particle &p2) {
        std::array<double, 3> add = {1., 0., 0.};
        p1.addToF(add);
        p2.addToF(add);
    });
    auto celllist = test->getCells();
    auto it1 = celllist[16].begin();
    auto it2 = celllist[11].begin();

    EXPECT_DOUBLE_EQ(it1->getF()[0], 6);
    EXPECT_DOUBLE_EQ(it2->getF()[0], 8);

}

TEST_F(LinkedCellTest, PeriodicTest1){
    test->addPeriodic(Boundary::TOP);
    test->addPeriodic(Boundary::BOTTOM);
    LennardJones lj{};
    std::shared_ptr<Container> test1 = test;
    lj.calculateF(test1);
    auto celllist = test->getCells();
    auto it1 = celllist[6].begin();
    EXPECT_DOUBLE_EQ( it1->getF()[0], -114.375);

    EXPECT_DOUBLE_EQ( it1->getF()[2], 0);

}

TEST_F(LinkedCellTest, PeriodicTest2){
    test->addPeriodic(Boundary::LEFT);
    test->addPeriodic(Boundary::RIGHT);
    LennardJones lj{};
    std::shared_ptr<Container> test1 = test;
    lj.calculateF(test1);
    auto celllist = test->getCells();
    auto it1 = celllist[6].begin();
    auto it2 = celllist[11].begin();
    EXPECT_DOUBLE_EQ( it1->getF()[1], -114.375);

    EXPECT_DOUBLE_EQ( it1->getF()[2], 0);
}

/**
 * @brief test applyF with several particles per cell
 */
TEST(LinkedCellTest_Outflow_Test, AppTest1){
    ParticleGenerator<LinkedCellContainer> cub{};
    std::shared_ptr<LinkedCellContainer> test1 = std::make_shared<LinkedCellContainer>();
    std::array<double, 3> domain{3., 3., 0.5};
    test1->setDomain(domain);
    test1->setSize(1., domain);
    cub.generateCuboid(test1, {.75, .75, 0}, {3, 3, 1}, .5, 1.0, {0., 1, 0.});

    test1->applyF([](Particle &p1, Particle &p2) {
        std::array<double, 3> add = {1., 0., 0.};
        p1.addToF(add);
        p2.addToF(add);
    });


    std::vector<ParticleContainer> celllist = test1->getCells();
    auto it1 = celllist[6].begin();
    auto it4 = celllist[12].begin();

    EXPECT_DOUBLE_EQ(it1->getF()[0], 8);
    EXPECT_DOUBLE_EQ(it4->getF()[0], 8);
}