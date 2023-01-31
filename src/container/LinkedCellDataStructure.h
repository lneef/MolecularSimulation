//
// Created by lukas on 11.01.23.
//
#pragma once



/**
 * @brief enum to store identifiers for sides
 */
enum class Boundary{LEFT, RIGHT, BOTTOM, TOP, FRONT, BACK};

#include "Reflecting.h"
#include "map"

/**
 * @brief upper class of all classes that implement the LinkedCell algorithm. It provides a uniform interface for the LinkedCell algorithm
 */
class LinkedCellDataStructure : public Container{
public:
     /**
     * @brief the given function to the particles in the container
     * @param fun function taking lvalue reference to particle
     */
    void apply(std::function<void(Particle &)> fun) override = 0;

     /**
     * @brief applies the given function to calculate the position of a particle
     * @param fun function taking lvalue reference to particle
     */
    void applyX(std::function<void(Particle &)> fun)  override = 0;

    /**
     * @brief overridden destructor to prevent memory leaks
     */
    ~LinkedCellDataStructure() override;

    /**
     * @brief calculates the force effective on particles using the given function
     * @param fun std::function taking two lvalue reference to particles
     */
    void applyF(std::function<void(Particle &, Particle &)> fun) override = 0;

    /**
     * @brief returns size of the container
     * @return number of particles stored in the container
     */
    size_t size() override = 0;

    /**
     * @brief adds a particle to linked cells
     * @param p rvalue reference to particle
     */
    void addParticle(Particle &&p) override = 0;

    /**
     * @brief function to add a periodic boundary
     * @param bound
     */
    static void addPeriodic(Boundary bound);

    /**
     * @brief set containing periodic boundaries
     */
    virtual void addParticle(Particle& p) = 0;
    /**
     * @brief removes all particles from the halo
     */

     /**
     * @brief funtion to initailize the linked cell algorithm for given domain and cutoff radius
     * @param rcutoff_arg cutoff radius for the linked cell algorithm
     * @param domain_arg domain the linked cells are covering
     * @param dim dimension of the simulation; currently only value 2 acceptable
     */
    virtual void setSize(double rcutoff_arg, std::array<double, 3> &domain_arg) = 0;

    /**
     * @brief function to obtain the domain
     * @return array of three doubles
     */
    virtual std::array<double, 3> &getDomain() = 0;

    /**
     * @brief remove all boundary conditions
     */
    static void clearBoundary();

    /**
     * @brief add reflecting boundary condition
     * @param bound boundary to which this condition applies
     * @param reflecting instance of Reflecting representing the boundary condition
     */
    static void addReflecting(Boundary bound, Reflecting &&reflecting);

protected:

    /**
     * @brief map containing reflecting boundary conditions
     */
    static std::map<Boundary,Reflecting> conditions;

    /**
     * @brief set containing periodic boundary conditions
     */
    static std::set<Boundary> periodic;

    /**
     * @brief cutoff radius of each cell in each dimension
     */
    static std::array<double, 3> cutoff;

    /**
     * @brief array to hold the number of cells in all three dimensions
     */
    static std::array<size_t, 3> mesh;

    /**
     * @brief vector spanning the domain which is covered by the cells
     */
    static std::array<double, 3> domain;

    /**
     * @brief implements c++20 contains for LinkedCellContainer for icpc
     * @param bound Boundary
     * @return true if for bound periodic boundary is specified, false otherwise
     */
    static bool containsPeriodic(Boundary bound);



};

