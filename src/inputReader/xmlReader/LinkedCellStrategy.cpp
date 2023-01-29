//
// Created by lukas on 11.01.23.
//

#include <iostream>
#include "LinkedCellStrategy.h"
#include "container/LinkedCellContainer.h"
#include "container/LinkedCellParallel.h"

namespace XMLReader {
    std::shared_ptr<LinkedCellDataStructure> &LinkedCellStrategy::get() {
        return lc;
    }

    std::shared_ptr<LinkedCellDataStructure> &LinkedCellStrategy::chose(size_t dim, std::string mode) {
        if(dim == 2){
            lc = std::make_shared<LinkedCellContainer>();
        }else{
            if(mode=="generic") {
                lc = std::make_shared<LinkedCell3D>();
            }else if(mode=="optimized"){
                lc = std::make_shared<LinkedCellParallel>();
            }
        }
        return lc;
    }
}

