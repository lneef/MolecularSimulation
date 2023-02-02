//
// Created by dominik on 12.12.22.
//

#pragma once

#include "../molsim-pskel.h"
#include "Simulation.h"
#include "utils//Thermostat.h"

namespace XMLReader {
    class temperature_pimpl : public XMLReader::temperature_pskel {
    private:
        /**
        * @brief Shared pointer of the simulation
        */
        std::shared_ptr<Simulation> sim;
        /**
        * @brief Current temperature
        */
        double temp;
        /**
         * @brief time step after which temperature values are updated
         */
        int n_thermo;
        /**
         * @brief Target temperature, which is the temperature you aim to reach
         */
        double temp_tar;
        /**
         * @brief Delta temperature, which is the maximal absolute temperature
         */
        double temp_del =  std::numeric_limits<double>::infinity();

        bool target_set = false;
    public:
        /**
         * @brief initializes the parser with a instance of Simulation
         * @param sim reference to shared pointer pointing to instance of Simulation
         */
        void init(std::shared_ptr<Simulation> &sim);

        /**
         * @brief Function that reads the initialization temperature
         */
        void temp_int(double) override;

        /**
         * @brief Function that reads the n_thermostat
         */
        void n_thermostat(int) override;

        /**
         * @brief Function that reads the target temperature
         */
        void temp_target(double) override;

        /**
         * @brief Function that reads the delta temperature
         */
        void temp_delta(double) override;

        /**
         * @brief Function that generates the thermostat pointer and passes the parameters to the simulation and the thermostat
         */
        void post_temperature() override;
    };
}