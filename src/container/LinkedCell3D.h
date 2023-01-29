//
// Created by lukas on 02.01.23.
//

#pragma once


#include <valarray>
#include "LinkedCellContainer.h"

class LinkedCell3D : public LinkedCellDataStructure {
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
    ~LinkedCell3D() override;

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
    LinkedCell3D();

    /**
     * @brief function to retrieve LinkedCellContainer at given index
     * @param i index of the LinkedCellContainer
     * @return reference to LinkedCellContainer
     */
    LinkedCellContainer &operator[](size_t i);

    void setSize(double cutOff_arg, std::array<double, 3> &domain_arg) override;

    std::array<double, 3> &getDomain() override;

    /**
     * @brief parallelized version of apply
     * @param fun function taking lvalue reference to Particle
     */
    void applyPar(std::function<void(Particle &)> fun) override;

    void setDomain(std::array<double, 3> &domain_arg);


private:
    std::vector<LinkedCellContainer> layers;

    std::array<double, 3> domain{};

    std::array<size_t, 3> mesh{};

    size_t layerSize = 0;

    size_t index(const std::array<double, 3> &pos) noexcept;

    void clearHalo();

    void preparePeriodic();

    void forceThreeD(Particle &p, size_t ind2D, size_t ind3D, std::function<void(Particle &, Particle &)> fun);

    void frontBackBoundary(double to_add, size_t ind, size_t oth);

    bool side(size_t ind3D);

    void update();

    void update(Particle &particle, size_t ind3D, size_t ind);

    size_t updatePeriodic(Particle &p, size_t ind3D);

    void applyFBoundary(std::function<void(Particle &, Particle &)> fun);

    void mirrorHorizontal(double dis, size_t i, LinkedCellContainer &counter);

    void mirrorVertical(double dis, size_t i, LinkedCellContainer &counter);

    void mirrorDiagonal(double dis, size_t i, LinkedCellContainer &counter);
};

