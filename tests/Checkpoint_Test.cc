#include "gtest/gtest.h"
#include <gmock/gmock-matchers.h>
#include "inputReader/xmlReader/XmlReader.h"

TEST(ReaderWriterTest, CheckpointTest) {
    std::string tes = "../tests/testinput/test2.xml";
    XMLReader::XmlReader xml{tes};
    std::shared_ptr<Simulation> sth1 = std::make_shared<Simulation>();
    std::shared_ptr<XMLReader::LinkedCellStrategy> lc1 = std::make_shared<XMLReader::LinkedCellStrategy>();
    xml.read(sth1, lc1);
    sth1->setParticle(lc1->get());
    sth1->checkpoint("checkpoint");

    std::string tes2 = "../tests/testinput/test_checkpoint.xml";
    XMLReader::XmlReader xml2{tes2};
    std::shared_ptr<Simulation> sth2 = std::make_shared<Simulation>();
    std::shared_ptr<XMLReader::LinkedCellStrategy> lc2 = std::make_shared<XMLReader::LinkedCellStrategy>();
    xml2.read(sth2, lc2);
    sth2->setParticle(lc2->get());

    std::vector<double> vel1;
    std::vector<double> vel2;

    lc1->get()->apply([&vel1](Particle& p) {
        std::array<double, 3> pos = p.getX();
        vel1.push_back(pos[0]);
        vel1.push_back(pos[1]);
        vel1.push_back(pos[2]);
    });


   lc2->get()->apply([&vel2](Particle& p) {
        std::array<double, 3> pos = p.getX();
        vel2.push_back(pos[0]);
        vel2.push_back(pos[1]);
        vel2.push_back(pos[2]);
    });

    for(int i = 0; i < vel1.size(); i++){
        EXPECT_THAT(vel1, testing::Pointwise(testing::DoubleEq(), vel2));
    }
}