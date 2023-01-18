//
// Created by lukas on 01.12.22.
//
#pragma once


#include "inputReader/xmlReader/molsim-pskel.h"
#include "simulation_pimpl.h"
#include "cuboid_input_pimpl.h"
#include "cuboid_pimpl.h"
#include "spheres_input_pimpl.h"
#include "boundaries_impl.h"
#include "temperature_pimpl.h"
#include "sphere_pimpl.h"
#include "membrane_pimpl.h"

namespace XMLReader {
    class molecular_pimpl : public molecular_pskel {
    public:

        void pre() override;

        /**
         * @brief function to process the subparser cuboid
         */
        void cuboid() override;

        /**
         * @brief function to process the subparser simulation
         */
        void simulation() override;

        /**
         * @brief function to process the subparser cuboid_input
         */
        void cuboid_input() override;

        /**
         * @brief function to process the subparser spheres_input
         */
        void spheres_input() override;

        /**
         * @brief function to process the subparser sphere
         */
        void sphere() override;

        /**
         * @brief function to process the subparser boundaries
         */
        void boundaries() override;

        /**
         * @brief function to process the subparser temperature
         */
        void temperature() override;

        /**
         * @brief function to process the subparser from_checkpoint
         */
        void from_checkpoint() override;

        void membrane() override;
        /**
         * @brief Function that generates the main parser molecular
         */
        void post_molecular() override;
    };

}