//
// Created by lukas on 28.01.23.
//

#include "BoundaryException.h"


const char *BoundaryException::what() const noexcept {
    return message.c_str();
}

BoundaryException::BoundaryException(std::string&& msg) {
    message = msg;

}
