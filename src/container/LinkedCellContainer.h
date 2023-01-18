#pragma once

#include "Container.h"
#include "ParticleContainer.h"
#include "Reflecting.h"
#include "LinkedCellDataStructure.h"
#include <vector>


class LinkedCell3D;
/**
 * @brief LinkedCellContainer implements the linked cell algorithm for a 2D simulation
 *
 * @warning The default boundary condition is outflow. If for a given boundary reflecting should be used, the condition must be added via addReflecting
 *
 * \image html linkedcell.png "Benchmark LinkedCellContainer" width=450cm
 * \image latex linkedcell.png "Benchmark LinkedCellContainer" width=10cm
 */

class LinkedCellContainer : public LinkedCellDataStructure {
public:


    /**
     * @brief the given function to the particles in the container
     * @param fun function taking lvalue reference to particle
     */
    void apply(std::function<void(Particle &)> fun) override;

    /**
     * @brief applies the given function to calculate the position of a particle
     * @param fun function taking lvalue reference to particle
     */
    void applyX(std::function<void(Particle &)> fun) override;

    /**
     * @brief overridden destructor to prevent memory leaks
     */
    ~LinkedCellContainer() override;

    /**
     * @brief calculates the force effective on particles using the given function
     * @param fun std::function taking two lvalue reference to particles
     */
    void applyF(std::function<void(Particle &, Particle &)> fun) override;

    /**
     * @brief returns size of the container
     * @return number of particles stored in the container
     */
    size_t size() override;

    /**
     * @brief adds a particle to linked cells
     * @param p rvalue reference to particle
     */
    void addParticle(Particle &&p) override;

    /**
     * @brief calculates the index of cells the given particle belongs to
     * @param p lvalue reference to particle
     * @return index of the particle in the linked cells data structure
     */
    static size_t index(Particle &p);

    /**
     * @brief default constructor for LinkedCellContainer
     */
    LinkedCellContainer();

    /**
     * @brief setter for cutoff radius
     * @param rcutoff_arg cutoff radius for the linked cell algorihtm
     */
    static void setRCutOff(double rcutoff_arg);

    /**
     * setter for the domain
     * @param domain_arg three dimensional vector representing the domain
     */
    static void setDomain(std::array<double, 3> &domain_arg);

    /**
     * @brief funtion to initailize the linked cell algorithm for given domain and cutoff radius
     * @param rcutoff_arg cutoff radius for the linked cell algorithm
     * @param domain_arg domain the linked cells are covering
     * @param dim dimension of the simulation; currently only value 2 acceptable
     */
     void setSize(double rcutoff_arg, std::array<double, 3> &domain_arg) override;


    /**
     * @brief function to get the cells for tests
     * @return std::vector representing the linked cells
     */
    std::vector<ParticleContainer>& getCells();

    /**
     * @brief returns the domain size
     * @return size of the domain
     */
    std::array<double, 3> &getDomain() override;

    /**
     * @brief returns particles in halo
     * @return lvalue reference to ParticleContainer containing particles in the halo
     */
    [[nodiscard]] const std::vector<std::reference_wrapper<ParticleContainer>> &getHalo() const;

    /**
     * @brief returns reference to boundary cells
     * @return std::vector containing reference wrappers to boundary cells
     */
    [[nodiscard]] const std::vector<std::reference_wrapper<ParticleContainer>> &getBoundary() const;

    /**
     * @brief set containing periodic boundaries
     */
    void addParticle(Particle& p) override;
    /**
     * @brief removes all particles from the halo
     */
    void clearHalo();

    static void setMesh(std::array<size_t, 3>& mesh_arg);

    ParticleContainer& operator[](size_t i);

    void clearBoundary();


private:

    /**
     * @brief std::vector to store the inner cells of the linked cell algorithm
     */
    std::vector<ParticleContainer> cells;

    /**
     * @brief array to hold the number of cells in all three dimensions
     */
    static std::array<size_t, 3> mesh;

    /**
     * @brief vector spanning the domain which is covered by the cells
     */
    static std::array<double, 3> domain;

    /**
     * @brief cutoff radius, i.e. length of the sides of the cells
     */
    static double rcutoff;
    /**
     * @brief updates the cell a particle is contained in after calculating the positions
     */
    void update();

    /**
     * @brief ParticleContainer storing the particles in the halo
     */
    std::vector<std::reference_wrapper<ParticleContainer>> halo;

    /**
     * @brief vector containing references to boundary cells
     */
    std::vector<std::reference_wrapper<ParticleContainer>> boundary;

    /**
     * @brief updates the boundary field after initialization
     */
    void setUpRef();

    void setUp();

    /**
     * @brief vector containing reelecting boundaries
     */


    /**
     * @brief applies reflecting boundary to particles in boundary cells
     * @param cond object representing reflecting boundary condition
     * @param fun function to calculate force
     */
    void applyFBoundary(Reflecting cond, std::function<void(Particle &, Particle &)> &fun);

    /**
     * @brief calculates the force between a particle and its right neighbours
     * @param i index of the cell of the particle
     * @param partial function which calculates the force between one specific particle and an arbitrary particle
     */
    void rightNeighbour(size_t i,const std::function<void(Particle &)> &partial);

    /**
     * @brief calculates the force between a particle and its upper neighbours
     * @param i index of the cell of the particle
     * @param partial function which calculates the force between one specific particle and an arbitrary particle
     */
    void upperNeighbour(size_t i, const std::function<void(Particle &)> &partial);

    /**
     * @brief calculates the force between a particle and its upper left neighbours
     * @param i index of the cell of the particle
     * @param partial function which calculates the force between one specific particle and an arbitrary particle
     */
    void upperLeftNeighbour(size_t i, const std::function<void(Particle &)> &partial);

     /**
     * @brief calculates the force between a particle and its upper right neighbours
     * @param i index of the cell of the particle
     * @param partial function which calculates the force between one specific particle and an arbitrary particle
     */
    void upperRightNeighbour(size_t i, const std::function<void(Particle &)> &partial);


    /**
     * @brief check if for given side periodic boundary is specified
     * @param ind index of a cell
     * @return true if periodic boundary is specified, false otherwise
     */
    bool side(size_t ind);

    /**
     * @brief mirrors particles in boundary cells in the halo if periodic boundary is specified
     * @param ind index of the cell the particle is contained in
     * @param p particle to be mirrored
     */
    void mirrorPeriodic(size_t ind, Particle& p);

    /**
     * @brief add a particle to the data structure without any checks
     * @param p rvalue reference to particle
     * @warning no checks on position are performed
     */
    void simpleAdd(Particle&& p);

    /**
     * @brief mirrors particle from halo cell with given index back into the inner cells
     * @param p lvalue reference to particle
     * @param ind index of the current halo cell
     * @return index of the cell particle was mirrored into
     */
    size_t mirror(Particle &p, size_t ind);

    /**
     * @brief update the cell a particle is contained in
     * @param p lvalue reference to particle
     * @param ind current index of the cell the particle is contained in
     */
    void update(Particle &p, size_t ind);


    /**
     * @brief updates the mirroring of particles int the boundary cells for periodic boundaty
     */
    void updatePeriodic();

    /**
     * @brief check test if cell with given index is in bottom boundary
     * @param ind index of a cell
     * @return true if cell in bottom boundary, false otherwise
     */
    static bool bottomBoundary(size_t ind);

    /**
     * @brief check test if cell with given index is in left boundary
     * @param ind index of a cell
     * @return true if cell in left boundary, false otherwise
     */
    static bool leftBoundary(size_t ind);

     /**
     * @brief check test if cell with given index is in right boundary
     * @param ind index of a cell
     * @return true if cell in right boundary, false otherwise
     */
    static bool rightBoundary(size_t ind);

     /**
     * @brief check test if cell with given index is in upper boundary
     * @param ind index of a cell
     * @return true if cell in upper boundary, false otherwise
     */
    static bool topBoundary(size_t ind);

    friend class LinkedCell3D;

    void forceTwoD(ParticleContainer &particles, size_t ind, std::function<void(Particle &, Particle &)> fun);


};
