//
// Created by dominik on 03.12.22.
//

#include "boundaries_impl.h"
#include "Reflecting.h"
#include "container/LinkedCellContainer.h"
#include <iostream>
#include "MolSimLogger.h"
#include "inputReader/xmlReader/LinkedCellStrategy.h"
#include "utils/BoundaryException.h"

namespace XMLReader {
    void boundaries_pimpl::top_boundary(const ::std::string &top) {
        MolSimLogger::logDebug("XMLReader: boundary {}", top);

        if (top == "reflecting"){
            LinkedCellDataStructure::addReflecting(Boundary::TOP, Reflecting(hor, cells->get()->getDomain()[1]));
        }else if(top == "periodic"){
            ++second;
            LinkedCellDataStructure::addPeriodic(Boundary::TOP);
        }
    }

    void boundaries_pimpl::bottom_boundary(const ::std::string &bottom) {
        MolSimLogger::logDebug("XMLReader: boundary {}", bottom);

        if (bottom == "reflecting"){
            LinkedCellDataStructure::addReflecting(Boundary::BOTTOM, Reflecting(hor, 0));
        }else if(bottom == "periodic"){
            ++second;
            LinkedCellDataStructure::addPeriodic(Boundary::BOTTOM);
        }
    }

    void boundaries_pimpl::left_boundary(const ::std::string &left) {
        MolSimLogger::logDebug("XMLReader: boundary {}", left);

        if (left == "reflecting"){
            LinkedCellDataStructure::addReflecting(Boundary::LEFT, Reflecting(vert, 0));
        }else if(left == "periodic"){
            ++first;
            LinkedCellDataStructure::addPeriodic(Boundary::LEFT);
        }
    }

    void boundaries_pimpl::right_boundary(const ::std::string &right) {
        MolSimLogger::logDebug("XMLReader: boundary {}", right);

        if (right == "reflecting"){
            LinkedCellDataStructure::addReflecting(Boundary::RIGHT, Reflecting(vert, cells->get()->getDomain()[0]));
        }else if(right == "periodic"){
            ++first;
            LinkedCellDataStructure::addPeriodic(Boundary::RIGHT);
        }
    }

    void boundaries_pimpl::post_boundaries() {
        std::stringstream stream;
        bool wrong_num= false;
        stream << "Periodic must be applied on bot boundaries in: ";
        if(first%2 == 1){
            stream<<"first";
            wrong_num = true;
        }
        if(second%2 == 1){
            stream<<"second";
            wrong_num = true;
        }
        if(third%2 == 1){
            stream<<"third";
            wrong_num = true;
        }

        if(wrong_num){
            stream<<"dimension"<<std::endl;
            throw BoundaryException(stream.str());
        }

    }

    void boundaries_pimpl::init(std::shared_ptr<LinkedCellStrategy> &cells_arg) {
        cells = cells_arg;
    }

    void boundaries_pimpl::front_boundary(const std::string & front) {
        if ( front == "reflecting"){
            LinkedCellDataStructure::addReflecting(Boundary::FRONT, Reflecting(dim3, 0));
        }else if(front == "periodic"){
            ++third;
            LinkedCellDataStructure::addPeriodic(Boundary::FRONT);
        }
    }

    void boundaries_pimpl::back_boundary(const std::string & back) {
        if (back == "reflecting"){
            LinkedCellDataStructure::addReflecting(Boundary::BACK, Reflecting(dim3, cells->get()->getDomain()[2]));
        }else if(back == "periodic"){
            ++third;
            LinkedCellDataStructure::addPeriodic(Boundary::BACK);
        }
    }
}