//
// Created by lukas on 11.01.23.
//
#pragma once


#include <memory>
#include "container/LinkedCellDataStructure.h"
#include "container/LinkedCell3D.h"

namespace XMLReader {
    class LinkedCellStrategy {
    public:
        std::shared_ptr<LinkedCellDataStructure> &get();

        std::shared_ptr<LinkedCellDataStructure> &chose(size_t dim, std::string mode);


    private:
        std::shared_ptr<LinkedCellDataStructure> lc;


    };
}