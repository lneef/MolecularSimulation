#pragma once
#include "limits"
#include "memory"
#include "container/Container.h"
/**
 * @brief Thermostat implements the thermostat functionality for the simulation
 * This includes cooling and heating the liquids
 */
class Thermostat {
private:
    /**
     * @brief Target temperature, which is the temperature you aim to reach
     */
    double temp_target;
    /**
     * @brief Delta temperature, which is the maximal absolute temperature
     */
    double temp_delta = std::numeric_limits<double>::infinity();
    /**
     * @brief Current temperature
     */
    double temp;

    /**
     * @brief calculates the scalar product for the velocity of the given Particles
     * @param p instance of Particle
     * @return scalar product of the velocity
     */
    static double scalarProduct(Particle& p);

public:
    /**
     * @brief applies the Thermostat to container of particles, i.e. a substance
     * @param particles shared pointer Container of particles
     */
    void applyThermostat(std::shared_ptr<Container>& particles);

    /**
     * @brief measure the temperature of a container of particles, i.e. a substance
     * @param particles shared pointer to Container of particles
     * @return
     */
    double measureTemp(std::shared_ptr<Container>& particles);

    /**
     * @brief constructor of the class Thermostat
     * @param current temperature
     */
    Thermostat(double current);

    /**
     * @brief getter for the field temp
     * @return double value equal to temp
     */
    double getTemp() const;

    /**
     * @brief setter for the field temp_target
     * @param target new target temperature
     */
    void setTarget(double target);

    /**
     * @brief setter for the field temp_delta
     * @param delta new max. temperature difference
     */
    void setDelta(double delta);
};

