//
// Created by dominik on 14.01.23.
//

#include "statistics_impl.h"

namespace XMLReader {
    void statistics_impl::init(std::shared_ptr<LinkedCellContainer> &lc, std::shared_ptr<Simulation> &simulation) {
        cells = lc;
        sim = simulation;
    }

    void statistics_impl::begin_rdf(int begin_rdf) {
        begin = begin_rdf;
    }

    void statistics_impl::end_rdf(int end_rdf) {
        end = end_rdf;
    }

    void statistics_impl::delta_rdf(double delta_rdf) {
        delta = delta_rdf;
    }

    void statistics_impl::n_statistics(int n_statistics) {
        n = n_statistics;
    }

    void statistics_impl::post_statistics() {
        std::shared_ptr<Statistics> statistics_p = std::make_shared<Statistics>(begin,end,delta);
        statistics_p->setParticles(cells);
        sim->setN_statistics(n);
        sim->setUse_statistics(true);
        sim->setStatistics(statistics_p);
    }
}