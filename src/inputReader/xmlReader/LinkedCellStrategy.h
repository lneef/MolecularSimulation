//
// Created by lukas on 11.01.23.
//
#pragma once


#include <memory>
#include "container/LinkedCellDataStructure.h"
#include "container/LinkedCell3D.h"

namespace XMLReader {

    /**
     * @brief class to determine which LinkedCellDatastructure is supposed to be used
     */
    class LinkedCellStrategy {
    public:

        /**
         * @brief getter for chosen data structure
         * @return lvalue reference to shared pointer to the chosen data structure
         */
        std::shared_ptr<LinkedCellDataStructure> &get();

        /**
         * @brief function that determines the desired LinkedCellDatastructure based on the input
         * @param dim dimension of the simulation
         * @param mode desired parallelization strategy
         * @return shared pointer to chosen LinkedCellDatastructure
         * @warning must be called before accessing any fields of the class
         */
        std::shared_ptr<LinkedCellDataStructure> &chose(size_t dim, int mode);


    private:

        /**
         * @brief chosen LinkedCellDatastructure is stored internally
         */
        std::shared_ptr<LinkedCellDataStructure> lc;


    };
}