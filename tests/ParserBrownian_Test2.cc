#include "gtest/gtest.h"
#include "inputReader/xmlReader/XmlReader.h"


/**
 * @brief test parser with temperature element and brownian motion
 */
TEST(ParserTestBrownian2, ParserTest) {
    std::string tes = "../tests/testinput/test3.xml";
    XMLReader::XmlReader xml{tes};
    std::shared_ptr<Simulation> sth = std::make_shared<Simulation>();
    std::shared_ptr<XMLReader::LinkedCellStrategy> lc = std::make_shared<XMLReader::LinkedCellStrategy>();
    xml.read(sth, lc);
    auto & force = sth->getForce();
    auto cont = dynamic_cast<LinkedCellContainer&>(*(lc->get()));
    auto &particles = cont.getCells();

    std::array<double, 3> velocity = particles[6].getParticles()[0].getV();
    EXPECT_NE(velocity[0], 0);
    EXPECT_NE(velocity[1], 0);
}