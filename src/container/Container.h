//
// Created by lukas on 29.11.22.
//

#pragma once

#include <functional>
#include "Particle.h"
/**
 * @brief Container provides a uniform interface for working with containers of particles
 */
class Container {
public:
    /**
     * @brief applies the given function to all elements of the container
     *
     * @param fun function to be applied to all elements
     */
    virtual void apply(std::function<void(Particle &)> fun) = 0;

    /**
     * @brief applies given function parallelized to all elements
     * @param fun function to be applied to all elements
     */
    virtual void applyPar(std::function<void(Particle&)> fun) = 0;

    /**
     * @brief calculates the position of the particles inside the container
     *
     * @param fun function to calculate the next position of a particles
     */
    virtual void applyX(std::function<void(Particle &)> fun) = 0;

    /**
     * @brief virtual destructor to prevent memory leaks
     */
    virtual ~Container();

    /**
     * @brief calculates the force acting upon a particle
     *
     * @param fun function taking two particles to calculate the force between them
     */
    virtual void applyF(std::function<void(Particle &, Particle &)> fun) = 0;

    /**
     * @brief calculates the number of particles stored in the container
     *
     * @return number of particles currently stored in the container
     */
    virtual size_t size() = 0;

    /**
     * @brief adds a particle to the conainer given as rvalue reference
     *
     * @param p rvalue reference to particle to be added to the container
     */
    virtual void addParticle(Particle&& p) = 0;


};


