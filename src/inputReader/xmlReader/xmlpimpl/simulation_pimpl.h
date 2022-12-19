//
// Created by lukas on 01.12.22.
//

#pragma once

#include <queue>
#include "../molsim-pskel.h"
#include "container/LinkedCellContainer.h"
#include "Simulation.h"

namespace XMLReader {
    class simulation_pimpl : public simulation_pskel {
    private:
        /**
         * @brief Container containing all linked cells for the simulation
         */
        std::shared_ptr<LinkedCellContainer> cells;
        /**
         * @brief Simulation to set the parameters of the simulation
         */
        std::shared_ptr<Simulation> sim;
        /**
         * @brief Domain size in each dimension
         */
        std::queue<double> domain{};
        /**
         * @brief Cutoff radius
         */
        double rCutOff;
    public:
        /**
         * @brief Function that initializes the container and the simulation
         */
        void init(std::shared_ptr<LinkedCellContainer> &cells_arg, std::shared_ptr<Simulation> &simu);
        /**
         * @brief Function that reads the end time of the simulation
         */
        void t_end(double) override;
        /**
         * @brief Function that reads the delta time of the simulation
         */
        void delta_t(double) override;
        /**
         * @brief Function that reads the domain size in x-dimension
         */
        void domain_size_x(double) override;
        /**
         * @brief Function that reads the domain size in y-dimension
         */
        void domain_size_y(double) override;
        /**
         * @brief Function that reads the domain size in z-dimension
         */
        void domain_size_z(double) override;
        /**
         * @brief Function that reads the cutoff radius
         */
        void cutOff_radius(double) override;
        /**
         * @brief Function that reads the name of the output name and sets it
         */
        void output_name(const ::std::string &) override;
        /**
         * @brief Function that reads the frequency of the output files and sets it
         */
        void output_frequency(int) override;

        /**
         * @brief changes the force calculation method of the simulation to LJGravitation
         */
        void g_gravitation(double ) override;
        /**
         * @brief Function that sets the domain of the cells
         */
        void post_simulation() override;
    };
}
