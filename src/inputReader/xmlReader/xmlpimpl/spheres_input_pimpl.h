//
// Created by dominik on 03.12.22.
//

#pragma once

#include "inputReader/xmlReader/molsim-pskel.h"
#include "../LinkedCellStrategy.h"
#include "../../Cuboid_file.h"

namespace XMLReader {
    class spheres_input_pimpl : public spheres_input_pskel {
    private:
        /**
         * @brief Container containing all linked cells for the simulation
         */
        std::shared_ptr<LinkedCellStrategy> cells;
        /**
         * @brief Path to the sphere txt-file
         */
        std::string path_spheres;
    public:
        /**
         * @brief Function that initializes the container
         */
        void init(std::shared_ptr<LinkedCellStrategy> &cell_arg);
        /**
         * @brief Function that reads the path of the file
         */
        void path(const ::std::string &) override;
        /**
         * @brief Function that reads the parameters of the spheres in the file and passes it to a file reader
         */
        void post_spheres_input() override;
    };
}