

#pragma once

#include<vector>
#include<iterator>
#include "Particle.h"
#include "Container.h"
#include "ParticleIterator.h"

/**
 * @brief Class ParticleContainer provides a container for storing particles.
 * 
 * It encapsulates the data structure used to store the particles.
 * To enable simple iteration over the stored particles, ParticleContainer implements the iterator pattern.
 *
 * \image html particlecontainer.png "Benchmark ParticleContainer" width=450cm
 * \image latex particlecontainer.png "Benchmark ParticleContainer" width=10cm
*/

class ParticleContainer : public Container {
private:
    /**
     * @brief The data structure where the particles are stored internally.
    */
    std::vector<Particle> particles;


public:
    /**
     * @brief Getter for field particles.
     *
     * @return vector of particles
    */
    std::vector<Particle> getParticles();


    void applyPar(std::function<void(Particle&)> fun) override;


    /**
     * @brief Setter for field particles
     *
     * @param particles vector of Particles
    */
    void setParticles(const std::vector<Particle> &particles);

    /**
     * @brief Inserts a Particle into ParticleContainer.
     *
     * @param particle rvalue reference to a Particle
    */
    void addParticle(Particle &&particle) override;

    /**
     * @brief Returns size of field particles.
     *
     * @return Number of particles stored in the respective ParticleContainer.
    */
    size_t size() override;

    /**
     * @brief overloaded access operator for ParticleContainer
     * @param i index which is accessed
     * @return lvalue reference to Particle at the passed index
     */
    Particle &operator[](size_t i);

    /**
     * @brief applies the given function to all particles and thus calculates the next positions
     *
     * @param fun function taking particles
     */
    void applyX(std::function<void(Particle &)> fun) override;

    /**
     * @brief applies the given function to all unordered pairs of particles
     *
     * @param fun function taking to particles as argument
     */
    void applyF(std::function<void(Particle &, Particle &)> fun) override;

    /**
     * @brief adds a particle to the container
     *
     * @param p lvalue reference to particle
     */
    void addParticle(Particle &p);

    /**
     * @bief applies the given function to all particles in the container
     *
     * @param fun function taking a particle that is applied to all particles
     */
    void apply(std::function<void(Particle &)> fun) override;

    /**
     * @brief override destructor to prevent memory leaks
     */
    ~ParticleContainer() override;

    /**
     * @brief deletes all particles in the container
     */
    void clear();

    /**
     * @brief returns a iterator pointing to the first element of the corresponding ParticleContainer
     * @return iterator to the first element
     */
    ParticleIterator begin();

    /**
     * @brief returns a iterator pointing behind the last element of the corresponding ParticleContainer
     * @return iterator pointing behind the last element
     */
    ParticleIterator end();

    /**
     * @brief removes the element pointed to by iterator
     * @param iterator iterator to element to be removed
     * @return iterator to element after the removed one
     */
    ParticleIterator remove(ParticleIterator &iterator);

};
