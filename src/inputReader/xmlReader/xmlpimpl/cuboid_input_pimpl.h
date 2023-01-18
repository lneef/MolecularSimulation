//
// Created by lukas on 01.12.22.
//
#pragma once

#include "inputReader/xmlReader/molsim-pskel.h"
#include "container/LinkedCellContainer.h"
#include "../../Cuboid_file.h"
#include "inputReader/xmlReader/LinkedCellStrategy.h"

namespace XMLReader {
    /**
     * @brief cuboid_input_pimpl implements the parser which reads a path to an cuboid txt-file
     */
    class cuboid_input_pimpl : public cuboid_input_pskel {
    private:
        /**
         * @brief Container containing all linked cells for the simulation
         */
        std::shared_ptr<LinkedCellStrategy> cells;
        /**
         * @brief Path to the cuboid txt-file
         */
        std::string path_cuboids;
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
         * @brief Function that reads the parameters of the cuboids in the file and passes it to a file reader
         */
        void post_cuboid_input() override;
    };

}

