#include "gtest/gtest.h"
#include "container/LinkedCell3D.h"
#include "inputReader/ParticleGenerator.h"


class LinkedCell3D_Test : public ::testing::Test{
protected:
    std::shared_ptr<LinkedCell3D> lc;
    void SetUp() override{
        ParticleGenerator<LinkedCell3D> gen{};
        lc= std::make_shared<LinkedCell3D>();
        std::array<double, 3> dom{3, 3, 3};
        lc->setSize(1, dom);
        gen.generateCuboid(lc, {0.5, 0.5, 0.5},{3, 3, 3}, 1, 1, {0, 0, 0});
    }
    void TearDown() override {
        lc->clearBoundary();
    }
};
TEST(LinkedCell3D_Test, AddTest){
    EXPECT_EQ(lc->size(), 27);
}