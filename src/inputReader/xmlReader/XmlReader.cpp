#include "XmlReader.h"
#include "MolSimLogger.h"

namespace XMLReader {
    XmlReader::XmlReader(std::string &s) {
        filename = s;
    }

    void XmlReader::read(std::shared_ptr<Simulation> &sim, std::shared_ptr<LinkedCellStrategy> &lc) {
        XMLReader::molecular_pimpl parser{};
        xml_schema::document doc_p(parser, "molecular");
        sim_p.parsers(double_p, double_p, double_p, double_p, double_p, double_p, string_p, int_p, double_p, double_p, int_p, string_p);
        temp_p.parsers(double_p, int_p, double_p, double_p);
        cub_in_p.parsers(string_p);
        cub_p.parsers(int_p, double_p, double_p, double_p, double_p, double_p, int_p, int_p, int_p, double_p, double_p,
                      double_p, double_p, double_p, bool_p);
        sph_in_p.parsers(string_p);
        sph_p.parsers(int_p, double_p, double_p, double_p, double_p, double_p, int_p, double_p, double_p, double_p,
                      double_p, double_p, bool_p);
        bou_in_p.parsers(string_p, string_p, string_p, string_p, string_p, string_p);
        fr_check.parsers(string_p);
        memb_p.parsers(double_p, double_p, double_p, double_p, double_p, int_p, int_p, int_p, double_p, double_p,
                       double_p, double_p, double_p, bool_p, double_p, double_p, double_p);
        stat_p.parsers(int_p,int_p,double_p,int_p);

        std::unique_ptr<Force> force = std::make_unique<LennardJones>();
        sim->setForce(force);

        sim_p.init(lc, sim);
        cub_p.init(lc, sim);
        cub_in_p.init(lc);
        sph_p.init(lc, sim);
        sph_in_p.init(lc);
        bou_in_p.init(lc);
        temp_p.init(sim);
        fr_check.init(lc);
        memb_p.init(lc, sim);
        stat_p.init(lc,sim);
        parser.parsers(cub_p, sim_p, cub_in_p, sph_p, sph_in_p, bou_in_p, temp_p, fr_check, memb_p,stat_p);
        parser.pre();
        try {
            MolSimLogger::logInfo("Parsing xml file:{}", filename);
            doc_p.parse(filename, xml_schema::flags::dont_validate);
        } catch (xml_schema::exception &e) {
            MolSimLogger::logError("{}", e.what());
            exit(-1);
        }
        parser.post_molecular();
    }

    void XmlReader::read(std::shared_ptr<Simulation> &sim){
        std::shared_ptr<LinkedCellStrategy> lc = std::make_shared<LinkedCellStrategy>();
        read(sim, lc);
    }

}
