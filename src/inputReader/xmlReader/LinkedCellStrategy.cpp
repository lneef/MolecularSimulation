//
// Created by lukas on 11.01.23.
//

#include <iostream>
#include "LinkedCellStrategy.h"
#include "container/LinkedCellContainer.h"

namespace XMLReader {
    std::shared_ptr<LinkedCellDataStructure> &LinkedCellStrategy::get() {
        return lc;
    }

    std::shared_ptr<LinkedCellDataStructure> &LinkedCellStrategy::chose(size_t dim) {
        if(dim == 2){
            lc = std::make_shared<LinkedCellContainer>();
        }else{
            lc = std::make_shared<LinkedCell3D>();
        }
        return lc;
    }
}

