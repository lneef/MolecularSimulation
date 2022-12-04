//
// Created by lukas on 01.12.22.
//
#pragma once

#include "inputReader/xmlReader/molsim-pskel.h"
#include "container/LinkedCellContainer.h"
#include "../../Cuboid_file.h"

namespace XMLReader {
    class cuboid_input_pimpl : public cuboid_input_pskel {
    private:
        std::shared_ptr<LinkedCellContainer> cells;
        std::string path_cuboids;
    public:
        void init(std::shared_ptr<LinkedCellContainer> &cell_arg);

        void path(const ::std::string &) override;

        void post_cuboid_input() override;
    };

}

