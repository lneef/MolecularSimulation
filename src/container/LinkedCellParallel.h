//
// Created by lukas on 29.01.23.
//

#pragma once


#include "LinkedCell3D.h"

class LinkedCellParallel : public LinkedCell3D{
public:
    /**
     * @brief applies the given function to all elements of the container
     *
     * @param fun function to be applied to all elements
     */
    void apply(std::function<void(Particle &)> fun) override;

    /**
     * @brief calculates the position of the particles inside the container
     *
     * @param fun function to calculate the next position of a particles
     */
    void applyX(std::function<void(Particle &)> fun) override;

    /**
     * @brief virtual destructor to prevent memory leaks
     */
    ~LinkedCellParallel() override;

    /**
     * @brief calculates the force acting upon a particle
     *
     * @param fun function taking two particles to calculate the force between them
     */
    void applyF(std::function<void(Particle &, Particle &)> fun) override;

    /**
     * @brief calculates the number of particles stored in the container
     *
     * @return number of particles currently stored in the container
     */
    size_t size() override;

    /**
     * @brief adds a particle given as rvalue reference to the container
     *
     * @param p rvalue reference to particle to be added to the container
     */
    void addParticle(Particle &&p) override;


    /**
     * @brief adds a particle given as rvalue reference to the container
     * @param p lvalue reference to particle to be added to the container
     */
    void addParticle(Particle &p) override;

    /**
     * @brief constructor
     */
    LinkedCellParallel();

    void setSize(double cutOff_arg, std::array<double, 3> &domain_arg) override;

    /**
     * @brief parallelized version of apply
     * @param fun function taking lvalue reference to Particle
     */
    void applyPar(std::function<void(Particle &)> fun) override;

};

