#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include "forceCalculation/Force.h"
#include "outputWriter/FileWriter.h"
#include "Simulation.h"


/**
 * @brief mocks the class Force
 */
class ForceMock : virtual public Force{
public:
    MOCK_METHOD(void , calculateF, (ParticleContainer &particles), (override));
};

/**
 * @brief mocks the class FileWriter
 */
class OutMock : public outputWriter::FileWriter{
public:
    MOCK_METHOD(void, plotParticles, (ParticleContainer &particles, const std::string &filename, int iteration) , (override));

};
/**
 * @brief test that calculation of f and output is called as often as expected
 */
TEST(SimulationTest, CallNumber){
    ParticleContainer par{};
    testing::Sequence s1, s2;
    std::unique_ptr<Force> force = std::make_unique<ForceMock>();
    std::unique_ptr<outputWriter::FileWriter> writer = std::make_unique<OutMock>();
    EXPECT_CALL(dynamic_cast<ForceMock&>(*force), calculateF).Times(12);
    EXPECT_CALL(dynamic_cast<OutMock&>(*writer), plotParticles).Times(1);
    Simulation sim(par, 1, 11, writer, force);
    sim.run();
}