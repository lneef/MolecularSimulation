//
// Created by lukas on 02.01.23.
//

#pragma once


#include <valarray>
#include "LinkedCellContainer.h"

/**
 * @brief LinkedCell3D implements the Linked Cell algorithm for three dimensions
 */
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
     * @brief constructor of LinkedCell3D
     */
    LinkedCell3D();

    /**
     * @brief function to retrieve LinkedCellContainer at given index
     * @param i index of the LinkedCellContainer
     * @return reference to LinkedCellContainer
     */
    LinkedCellContainer &operator[](size_t i);


    /**
     * @brief function to initialize the data structure
     * @param cutOff_arg cutoff radius
     * @param domain_arg domain of the simulation
     */
    void setSize(double cutOff_arg, std::array<double, 3> &domain_arg) override;

    std::array<double, 3> &getDomain() override;

    /**
     * @brief parallelized version of apply
     * @param fun function taking lvalue reference to Particle
     */
    void applyPar(std::function<void(Particle &)> fun) override;

    /**
     * @brief setter for the domain
     * @param domain_arg array of 3 double
     */
    void setDomain(std::array<double, 3> &domain_arg);


protected:
    /**
     * @brief vector containing the layers of the three dimensional Grid
     */
    std::vector<LinkedCellContainer> layers;

    /**
     * @brief number of cells of each two dimensional grid
     */
    size_t layerSize = 0;

    /**
     * @brief calculates the index of the layer of particle belongs to
     * @param pos current position of the particle
     * @return index of the layer
     */
    static size_t index(const std::array<double, 3> &pos) noexcept;

    /**
     * @brief remove all particles in the halo cells
     */
    void clearHalo();

    /**
     * @brief function to mirror particles for periodic boundary condition
     */
    void preparePeriodic();

    /**
     * @brief calculates force between a given particle and all neighbouring cells in the next layer
     * @param p particle that is considered
     * @param ind2D index of the particle in its current layer
     * @param ind3D index of the layer the particle belongs to
     * @param fun force calculation routine
     */
    void forceThreeD(Particle &p, size_t ind2D, size_t ind3D, std::function<void(Particle &, Particle &)> fun);

    /**
     * @brief mirror the particles in the third dimension
     * @param to_add distance to the cell on the opposite side
     * @param ind index of the layer a particle belongs to
     * @param oth index of the halo layer on the opposite side
     */
    void frontBackBoundary(double to_add, size_t ind, size_t oth);

    /**
     * @brief determine if periodic boundary is applicable for given layer
     * @param ind3D index of a layer
     * @return
     */
    bool side(size_t ind3D);

    /**
     * @brief update the cell a particle is contained in after each iteration
     */
    void update();

    /**
     * @brief updated the indices of the given particle and add it to new location
     * @param particle particle that is considered
     * @param ind3D index of the layer the particle is contained in
     * @param ind index of the particle in the layer
     */
    void update(Particle &particle, size_t ind3D, size_t ind);

    /**
     * @brief move particle at a periodic boundary
     * @param p particle that is considered
     * @param ind3D index of the layer a particle is contained in
     * @return new index of the particle in the third dimension
     */
    size_t updatePeriodic(Particle &p, size_t ind3D);

    /**
     * @brief apply reflecting boundary condition in third dimension
     * @param fun force calculation routine
     */
    void applyFBoundary(std::function<void(Particle &, Particle &)> fun);

    /**
     * @brief mirror particles in the inner cells
     * @param dis length of the domain in the third dimension
     * @param i index of the layer a particle belongs to
     * @param counter layer on the opposite side
     */
    void mirrorHorizontal(double dis, size_t i, LinkedCellContainer &counter);

    /**
     * @brief mirror particle in the boundary cells
     * @param dis length of the domain in the third dimension
     * @param i index of the layer a particle belongs to
     * @param counter layer on the opposite side
     */
    void mirrorVertical(double dis, size_t i, LinkedCellContainer &counter);

    /**
     * @brief mirror particle in the corner cells
     * @param dis length of the domain in the third dimension
     * @param i index of the layer a particle belongs to
     * @param counter layer on the opposite side
     */
    void mirrorDiagonal(double dis, size_t i, LinkedCellContainer &counter);
};

