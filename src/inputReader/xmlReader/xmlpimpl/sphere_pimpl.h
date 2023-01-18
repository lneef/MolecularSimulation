//
// Created by dominik on 14.12.22.
//

#pragma once

#include <queue>
#include "inputReader/xmlReader/molsim-pskel.h"
#include "Simulation.h"
#include "inputReader/xmlReader/LinkedCellStrategy.h"

namespace XMLReader {
    class sphere_pimpl : public sphere_pskel {
    private:
        /**
         * @brief shared pointer to LinkedCellContainer where particles are stores
         */
        std::shared_ptr<LinkedCellStrategy> cells;

        /**
         * @brief simulation that works with the sphere
         */
        std::shared_ptr<Simulation> sim;

        /**
         * @brief center of the sphere 2D
         */
        std::queue<double> centre;

        /**
         * @brief radius of the sphere
         */
        int radius_r;

        /**
         * @brief mesh width of the sphere
         */
        double width;

        /**
         * @brief mass of each particle
         */
        double mass_r;

        /**
         * @brief velocity of each particle in all 2 dimensions
         */
        std::queue<double> vel;

        /**
         * @brief flag for initialization with brownian motion
         */
        bool browMot = true;

        /**
         * @brief type of the particle
         */
        int type_s = 1;

        /**
         * @brief sigma value of the particles
         */
        double sigma_s = 1;

        /**
         * @brief epsilon value of the particles
         */
        double epsilon_s = 5;
    public:

        /**
         * @brief initializes the parser
         * @param lc reference to shared pointer pointing to a LinkedCellContainer
         * @param sim reference to shared pointer pointing to instance of Simulation
         */
        void init(std::shared_ptr<LinkedCellStrategy> &lc, std::shared_ptr<Simulation> &sim);

        /**
         * @brief function to process x coordinate of the center
         */
        void center_x(double) override;

        /**
         * @brief function to process y coordinate of the center
         */
        void center_y(double) override;

        /**
         * @brief function to process z coordinate of the center
         */
        void center_z(double) override;

        /**
         * @brief function to process the radius
         */
        void radius(int) override;

        /**
         * @brief function to process the mesh width
         */
        void mesh_width(double) override;

        /**
         * @brief function to process the mass
         */
        void mass(double) override;

        /**
         * @brief function to process the velocity in x direction
         */
        void velocity_x(double) override;

        /**
         * @brief function to process the velocity in y direction
         */
        void velocity_y(double) override;

        /**
         * @brief function to process the velocity in z direction
         */
        void velocity_z(double) override;

        /**
         * @brief function to process the flag for brownian motion
         */
        void brownianMotion(bool) override;

        /**
         * @brief function to process the type of the particles
         */
        void type(int) override;

        /**
         * @brief function to process the sigma value(Lennard Jones) for the particles
         */
        void sigma(double) override;

        /**
         * @brief function to process the epsilon value(Lennard Jones) for the particles
         */
        void epsilon(double) override;

        /**
         * @brief generates sphere and stores particles in cells
         */
        void post_sphere() override;
    };
}
