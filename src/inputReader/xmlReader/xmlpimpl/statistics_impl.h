//
// Created by dominik on 14.01.23.
//

#pragma once

#include "../molsim-pskel.h"
#include "Simulation.h"
#include "inputReader/xmlReader/LinkedCellStrategy.h"
#include "../../../Statistics.h"

namespace XMLReader{
     class statistics_impl : public XMLReader::statistics_pskel{
     private:
         /**
          * @brief First value of the rdf distance intervals
          */
         int begin;
         /**
          * @brief Last value of the rdf distance intervals
          */
         int end;
         /**
          * @brief Delta r, which defines the length of the rdf distance intervals
          */
         double delta;
         /**
          * @brief Time step after which statistics will be calculated
          */
         int n;
         /**
          * @brief Shared pointer of the simulation
          */
         std::shared_ptr<Simulation> sim;
         /**
          * @brief Shared pointer of the particles in the simulation
          */
         std::shared_ptr<LinkedCellStrategy> cells;
     public:
         /**
         * @brief initializes the parser
         * @param Reference to shared pointer pointing to the particles
         * @param Reference to shared pointer pointing to the simulation
         */
         void init(std::shared_ptr<LinkedCellStrategy> &lc,std::shared_ptr<Simulation> &sim);
         /**
         * @brief Function that reads the first value of the rdf distance intervals
         */
         void begin_rdf (int ) override;
         /**
         * @brief Function that reads the last value of the rdf distance intervals
         */
         void end_rdf (int ) override;
         /**
         * @brief Function that reads the delta r
         */
         void delta_rdf (double ) override;
         /**
         * @brief Function that reads the time step after which statistics will be calculated
         */
         void n_statistics (int ) override;
         /**
         * @brief Function that generates the statistics pointer and passes the parameters to the simulation and the statistics
         */
         void post_statistics() override;
     };
}