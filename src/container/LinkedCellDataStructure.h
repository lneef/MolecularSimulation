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
     * @brief adds a reflecting boundary condition to the container
     * @param ref rvalue reference to object of type Reflecting
     */
    static void addReflecting(Reflecting &&ref);

    static void addThird(Boundary bound, Reflecting&& ref);

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

    virtual std::array<double, 3> &getDomain() = 0;

    static void clearBoundary();

protected:
    static std::vector<Reflecting> conditions;

    static std::map<Boundary,Reflecting> frontBack;

    static std::set<Boundary> periodic;

    /**
     * @brief implements c++20 contains for LinkedCellContainer for icpc
     * @param bound Boundary
     * @return true if for bound periodic boundary is specified, false otherwise
     */
    static bool containsPeriodic(Boundary bound);



};

