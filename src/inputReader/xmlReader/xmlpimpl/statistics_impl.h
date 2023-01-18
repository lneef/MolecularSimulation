//
// Created by dominik on 14.01.23.
//

#pragma once

#include "../molsim-pskel.h"
#include "Simulation.h"
#include "Statistics.h"
#include "container/LinkedCellContainer.h"

namespace XMLReader{
     class statistics_impl : public XMLReader::statistics_pskel{
     private:
         int begin;
         int end;
         double delta;
         int n;
         std::shared_ptr<Simulation> sim;
         std::shared_ptr<LinkedCellContainer> cells;
     public:
         void init(std::shared_ptr<LinkedCellContainer> &lc,std::shared_ptr<Simulation> &sim);
         void begin_rdf (int ) override;
         void end_rdf (int ) override;
         void delta_rdf (double ) override;
         void n_statistics (int ) override;
         void post_statistics() override;
     };
}