#pragma once
#include "Force.h"
#include "container/ParticleContainer.h"

/**
 * @brief SLennardJones implements the interface provided by Force
 *
 * SLennardJones implements the interface provided by Force for calculating the force between SLennardJones-particles
*/
class SLennardJones: public Force {
public:

    /**
     * @brief calculates the force effective on the particles inside the given container
     * @param particles Container holding the particles for which the force is calculated
     */
    void calculateF(std::shared_ptr<Container>& particles) override;

    /**
     * @brief calculates the force between two particles based on the harmonic potential
     * @param p1 first particle of the interaction
     * @param p2 second particle of the interaction
     */
    static void calculateF(Particle& p1, Particle& p2, double rc, double rl);

    ~SLennardJones() override;

    /**
     * @brief constructor of SLennardJones
     */
     SLennardJones(double rc_arg, double rl_arg);

private:

    /**
     * @brief the cutoff radius
     */
    double rc;

    /**
     * @brief radius-low
     */
    double rl;
};
