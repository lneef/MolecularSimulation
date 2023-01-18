//
// Created by lukas on 01.12.22.
//
#pragma once

#include <queue>
#include "inputReader/xmlReader/molsim-pskel.h"
#include "inputReader/xmlReader/LinkedCellStrategy.h"
#include "Simulation.h"
#include "forceCalculation/MembraneForce.h"

namespace XMLReader {
    class membrane_pimpl : public membrane_pskel {
    private:
        /**
         * @brief Container containing all linked cells for the simulation
         */
        std::shared_ptr<LinkedCellStrategy> cells;
        /**
         * @brief Queue containing the position of the cuboid
         */
        std::queue<double> pos;
        /**
         * @brief Queue containing the quantity of molecules of the cuboid
         */
        std::queue<int> num;
        /**
         * @brief Queue containing the velocity of the cuboid
         */
        std::queue<double> vel;
        /**
         * @brief Mass of the molecules in the cuboid
         */
        double m;
        /**
         * @brief Distance between the molecules in the cuboid
         */
        double width;

        /**
         * @brief simulation which works with the particles
         */
        std::shared_ptr<Simulation> sim;

        /**
         * @brief flag for initialization with brownian motion, default: true to enable backwards compatability
         */
        bool browMot = true;

        /**
         * @brief sigma value of the particles for Lennard Jones potential
         */
        double sigma_p = 1;

        /**
         * @brief epsilon value of the particles for Lennard Jones potential
         */
        double epsilon_p = 5;

        double f = 1;

        double stiff_const;

        double bond_len;
    public:
        /**
         * @brief Function that initializes the container
         */
        void init(std::shared_ptr<LinkedCellStrategy> &lc, std::shared_ptr<Simulation> &sim_arg);
        /**
         * @brief Function that reads the position in x-dimension
         */
        void lower_left_x(double) override;
        /**
         * @brief Function that reads the position in y-dimension
         */
        void lower_left_y(double) override;
        /**
         * @brief Function that reads the position in z-dimension
         */
        void lower_left_z(double) override;
        /**
         * @brief Function that reads the quantity in x-dimension
         */
        void number_x(int) override;
        /**
         * @brief Function that reads the quantity in x-dimension
         */
        void number_y(int) override;
        /**
         * @brief Function that reads the quantity in x-dimension
         */
        void number_z(int) override;
        /**
         * @brief Function that reads the distance between the molecules
         */
        void mesh_width(double) override;
        /**
         * @brief Function that reads mass of the molecules
         */
        void mass(double) override;
        /**
         * @brief Function that reads the velocity in x-dimension
         */
        void velocity_x(double) override;
        /**
         * @brief Function that reads the velocity in x-dimension
         */
        void velocity_y(double) override;
        /**
         * @brief Function that reads the velocity in x-dimension
         */
        void velocity_z(double) override;

        /**
         * @brief function to process to flag for brownian motion
         */
        void brownianMotion(bool) override;

        /**
         * @brief function to process the sigma value for the particle
         */
        void sigma(double ) override;

        /**
         * @brief function to process the epsilon value for the particles
         */
        void epsilon(double ) override;

        void fz_up(double ) override;

        void stiffness_const(double ) override;

        void bond_length(double ) override;

        /**
         * @brief Function that generates the cuboids
         */
        void post_membrane() override;
    };

}

