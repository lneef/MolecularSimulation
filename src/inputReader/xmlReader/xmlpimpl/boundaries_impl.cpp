//
// Created by dominik on 03.12.22.
//

#include "boundaries_impl.h"
#include "Reflecting.h"
#include "container/LinkedCellContainer.h"
#include <iostream>
#include "MolSimLogger.h"
#include "inputReader/xmlReader/LinkedCellStrategy.h"

namespace XMLReader {
    void boundaries_pimpl::top_boundary(const ::std::string &top) {
        MolSimLogger::logDebug("XMLReader: boundary {}", top);

        if (top == "reflecting"){
            LinkedCellDataStructure::addReflecting(Reflecting(hor, cells->get()->getDomain()[1]));
        }else if(top == "periodic"){
            LinkedCellDataStructure::addPeriodic(Boundary::TOP);
        }
    }

    void boundaries_pimpl::bottom_boundary(const ::std::string &bottom) {
        MolSimLogger::logDebug("XMLReader: boundary {}", bottom);

        if (bottom == "reflecting"){
            LinkedCellDataStructure::addReflecting(Reflecting(hor, 0));
        }else if(bottom == "periodic"){
            LinkedCellDataStructure::addPeriodic(Boundary::BOTTOM);
        }
    }

    void boundaries_pimpl::left_boundary(const ::std::string &left) {
        MolSimLogger::logDebug("XMLReader: boundary {}", left);

        if (left == "reflecting"){
            LinkedCellDataStructure::addReflecting(Reflecting(vert, 0));
        }else if(left == "periodic"){
            LinkedCellDataStructure::addPeriodic(Boundary::LEFT);
        }
    }

    void boundaries_pimpl::right_boundary(const ::std::string &right) {
        MolSimLogger::logDebug("XMLReader: boundary {}", right);

        if (right == "reflecting"){
            LinkedCellDataStructure::addReflecting(Reflecting(vert, cells->get()->getDomain()[0]));
        } if(right == "periodic"){
            LinkedCellDataStructure::addPeriodic(Boundary::RIGHT);
        }
    }

    void boundaries_pimpl::post_boundaries() {
    }

    void boundaries_pimpl::init(std::shared_ptr<LinkedCellStrategy> &cells_arg) {

    }

    void boundaries_pimpl::front_boundary(const std::string & front) {
        if ( front == "reflecting"){
            LinkedCellDataStructure::addReflecting(Reflecting(vert, 0));
        }else if(front == "periodic"){
            LinkedCellDataStructure::addPeriodic(Boundary::FRONT);
        }
    }

    void boundaries_pimpl::back_boundary(const std::string & back) {
        if (back == "reflecting"){
            LinkedCellDataStructure::addReflecting(Reflecting(vert, cells->get()->getDomain()[2]));
        }else if(back == "periodic"){
            LinkedCellDataStructure::addPeriodic(Boundary::BACK);
        }
    }
}