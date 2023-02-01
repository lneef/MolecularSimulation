#pragma once
#include "Force.h"
#include "container/ParticleContainer.h"

/**
 * @brief MembraneForce implements the interface provided by Force
 *
 * MembraneForce implements the interface provided by Force for calculating the force between membrane-particles
*/
class MembraneForce: public Force {
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
    static void calculateF(Particle& p1, Particle& p2, double r0, double k);

    ~MembraneForce() override;

    static void membraneForce(std::shared_ptr<Container>& particles, double f_up, double g, std::unique_ptr<Force>& force, double time);

    /**
     * @brief constructor of MembraneForce
     */
     MembraneForce(double r0_arg, double k_arg);

private:

    /**
     * @brief average bond length of a molecule pair
     */
    double r0;

    /**
     * @brief stiffness constant
     */
    double k;
};
