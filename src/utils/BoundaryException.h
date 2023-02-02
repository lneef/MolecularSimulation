//
// Created by lukas on 28.01.23.
//
#pragma once


#include <exception>
#include <string>

/**
 * @brief Execption that is thrown is periodic boundries are specified incorrectly
 */
class BoundaryException : public std::exception{
public:

    /**
     * @brief overriden what function
     */
    [[nodiscard]] const char*
    what() const noexcept override;

    /**
    * @brief constructor of BoundaryException
    * @param msg rvalue reference to string thar indicates the error
    */
    BoundaryException(std::string&& msg);

private:

    /**
     * @brief error message
     */
    std::string message;

};


