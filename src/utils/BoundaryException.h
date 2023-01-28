//
// Created by lukas on 28.01.23.
//
#pragma once


#include <exception>
#include <string>

class BoundaryException : public std::exception{
public:
    [[nodiscard]] const char*
    what() const noexcept override;

    BoundaryException(std::string&& msg);

private:
    std::string message;

};


