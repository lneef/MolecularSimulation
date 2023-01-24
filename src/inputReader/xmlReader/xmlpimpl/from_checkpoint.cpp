//
// Created by lukas on 18.12.22.
//

#include "from_checkpoint.h"
#include "../../CheckPointReader.h"
#include "inputReader/xmlReader/LinkedCellStrategy.h"

namespace XMLReader {
    void from_checkpoint_pimpl::init(std::shared_ptr<LinkedCellStrategy> &cells_arg) {
        this->cells = cells_arg;

    }

    void from_checkpoint_pimpl::path(const std::string &s) {
        filename = s;
    }

    void from_checkpoint_pimpl::post_from_checkpoint() {
        inputReader::CheckPointReader checkpoint_input{filename};
        checkpoint_input.readInto(cells->get());

    }
}