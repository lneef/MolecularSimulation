#pragma once

#include <memory>
#include <iostream>
#include "../utils/MaxwellBoltzmannDistribution.h"
#include "container/ParticleContainer.h"
#include "utils/ArrayUtils.h"

/**
 * @brief ParticleGenerator generates cuboids of particles form input values
 *
 * @param T template argument for container type that is used
 *
 */
template<typename T>
class ParticleGenerator {
private:

    /**
     * @brief dimension for the brownian motion
     */
    size_t dim = 2;

    /**
     * @brief default mean for brownian motion
     */
    double meanV = 0.1;

    /**
     * @brief mb is used to obtain random values according to the MaxwellBoltzmann distribution
     */
    MaxwellBoltzmannDistribution mb;

    /**
     * @brief Calculates velocity of particles with given velocity and the MaxwellBoltzmannDistribution
     *
     * @param v Array representation of the velocities of the particles in all three dimensions
     * @return velocity with thermal friction
    */
    std::array<double, 3> calculateV(std::array<double, 3> v) {
        //calculate velocity with thermal friction
        std::array<double, 3> mbV = mb.maxwellBoltzmannDistributedVelocity(meanV, dim);
        std::array<double, 3> newV{};
        for (int i = 0; i < 3; i++) {
            newV[i] = v[i] + mbV[i];
        }
        return newV;
    }

    /**
     * @brief Calculates velocity of particles with given mean of the velocity and the MaxwellBoltzmannDistribution
     *
     * @param meanVelocity Double value representation of mean of the velocity
     * @return velocity with thermal friction
    */
    std::array<double, 3> calculateVWithMean(double meanVelocity) {
        //calculate velocity with thermal friction
        std::array<double, 3> mbV = mb.maxwellBoltzmannDistributedVelocity(meanVelocity, 2);
        return mbV;
    }

public:
    explicit ParticleGenerator(size_t dim_arg = 2){
        dim = dim_arg;
    }

    virtual ~ParticleGenerator() = default;

    /**
     * @brief Generates cuboids of particles with given values and stores them in instance of ParticleContainer
     *
     * @param particles ParticleContainer where the particles are stored
     * @param x Array representation three dimensional position of the particle
     * @param n Array representation of quantity of particles in each dimension
     * @param h Double representation of distance between the particles
     * @param m Double representation of mass of one particle
     * @param v Array representation of the velocities of the particles in all three dimensions
     * @param sigma sigma value for Lennard Jones
     * @param epsilon epsilon value for Lennard Jones
     * @param type type of the particle
    */
    void
    generateCuboid(std::shared_ptr<T> &particles, std::array<double, 3> x, std::array<int, 3> n, double h, double m,

                   std::array<double, 3> v, double sigma = 1, double epsilon = 5, int type=1) {
        std::array<double, 3> newX{};
        //iterate over cuboid in each dimension
        for (int x_cord = 0; x_cord < n[0]; x_cord++) {
            newX[0] = x[0] + (x_cord * h);
            for (int y_cord = 0; y_cord < n[1]; y_cord++) {
                newX[1] = x[1] + (y_cord * h);
                for (int z_cord = 0; z_cord < n[2]; z_cord++) {
                    newX[2] = x[2] + (z_cord * h);
                    std::array<double, 3> newV = calculateV(v);
                    //stores particle in ParticleContainer
                    particles->addParticle(Particle(newX, newV, m, sigma, epsilon, type));
                }
            }
        }
    }

    /**
     * @brief Generates cuboids of particles with given values and the MaxwellBoltzmannDistribution and stores them in instance of ParticleContainer
     * @param particles ParticleContainer where the particles are stored
     * @param x Array representation three dimensional position of the particle
     * @param n Array representation of quantity of particles in each dimension
     * @param v Array representation of the velocities of the particles in all three dimensions
     * @param h Double representation of distance between the particles
     * @param m Double representation of mass of one particle
     * @param meanVelocity Double representation of mean of the velocity
     * @param sigma sigma value for Lennard Jones
     * @param epsilon epsilon value for Lennard Jones
     * @param type type of the particle
     */
    void
    generateCuboidBrownian(std::shared_ptr<T> &particles, std::array<double, 3> x, std::array<int, 3> n,
                           std::array<double, 3> v, double h, double m,
                           double meanVelocity, double sigma = 1, double epsilon = 5, int type=1) {
        std::array<double, 3> newX{};
        //iterate over cuboid in each dimension
        for (int x_cord = 0; x_cord < n[0]; x_cord++) {
            newX[0] = x[0] + (x_cord * h);
            for (int y_cord = 0; y_cord < n[1]; y_cord++) {
                newX[1] = x[1] + (y_cord * h);
                for (int z_cord = 0; z_cord < n[2]; z_cord++) {
                    newX[2] = x[2] + (z_cord * h);
                    std::array<double, 3> mbV = calculateVWithMean(meanVelocity);
                    std::array<double, 3> newV{};
                    for (int i = 0; i < 3; i++) {
                        newV[i] = v[i] + mbV[i];
                    }
                    //stores particle in ParticleContainer
                    particles->addParticle(Particle(newX, newV, m, sigma, epsilon, type));
                }
            }
        }
    }

    /**
     * @brief Generates mebrane of particles with given values and the MaxwellBoltzmannDistribution and stores them in instance of ParticleContainer
     * @param particles ParticleContainer where the particles are stored
     * @param x Array representation three dimensional position of the particle
     * @param n Array representation of quantity of particles in each dimension
     * @param v Array representation of the velocities of the particles in all three dimensions
     * @param h Double representation of distance between the particles
     * @param m Double representation of mass of one particle
     * @param meanVelocity Double representation of mean of the velocity
     * @param sigma sigma value for Lennard Jones
     * @param epsilon epsilon value for Lennard Jones
     * @param type type of the particle
     * @param fz force in z up
     */
    void
    generateMembraneBrownian(std::shared_ptr<T> &particles, std::array<double, 3> x, std::array<int, 3> n,
                           std::array<double, 3> v, double h, double m,
                           double meanVelocity, double sigma = 1, double epsilon = 5, int type=3, double fz = 1) {
        std::array<double, 3> newX{};
        //iterate over cuboid in each dimension
        std::array<int, 2> index;
        for (int x_cord = 0; x_cord < n[0]; x_cord++) {
            newX[0] = x[0] + (x_cord * h);
            for (int y_cord = 0; y_cord < n[1]; y_cord++) {
                newX[1] = x[1] + (y_cord * h);
                for (int z_cord = 0; z_cord < n[2]; z_cord++) {
                    newX[2] = x[2] + (z_cord * h);
                    std::array<double, 3> mbV = calculateVWithMean(meanVelocity);
                    std::array<double, 3> newV{};
                    for (int i = 0; i < 3; i++) {
                        newV[i] = v[i] + mbV[i];
                    }
                    //stores particle in ParticleContainer
                    Particle p = Particle(newX, newV, m, sigma, epsilon, type);
                    index[0] = x_cord;
                    index[1] = y_cord;
                    p.setIndex(index);
                    std::array<double, 3> force = {0,0,fz};
                    if((index[0] == 17 && index[1] == 24) ||( index[0] == 17 && index[1] == 25 )|| (index[0] == 18 && index[1] == 24 )||( index[0] == 18 && index[1] == 25)){
                        p.setF(force);
                    }
                    particles->addParticle(p);
                }
            }
        }
    }

    /**
     * @brief Generates cuboids of particles with given values and stores them in instance of ParticleContainer
     * @param particles ParticleContainer where the particles are stored
     * @param x Array representation three dimensional position of the particle
     * @param n Array representation of quantity of particles in each dimension
     * @param v Array representation of the velocities of the particles in all three dimensions
     * @param h Representation of distance between the particles
     * @param m Representation of mass of one particle
     * @param sigma sigma value for Lennard Jones
     * @param epsilon epsilon value for Lennard Jones
     * @param type type of the particle
     */
    void
    generateCuboidNoBrownian(std::shared_ptr<T> &particles, std::array<double, 3> x, std::array<int, 3> n,
                             std::array<double, 3> v, double h, double m, double sigma = 1, double epsilon = 5, int type= 1) {
        std::array<double, 3> newX{};
        //iterate over cuboid in each dimension
        for (int x_cord = 0; x_cord < n[0]; x_cord++) {
            newX[0] = x[0] + (x_cord * h);
            for (int y_cord = 0; y_cord < n[1]; y_cord++) {
                newX[1] = x[1] + (y_cord * h);
                for (int z_cord = 0; z_cord < n[2]; z_cord++) {
                    newX[2] = x[2] + (z_cord * h);
                    //stores particle in ParticleContainer
                    particles->addParticle(Particle(newX, v, m, sigma, epsilon, type));
                }
            }
        }
    }

    /**
     * @brief Generates membrane of particles with given values and stores them in instance of ParticleContainer
     * @param particles ParticleContainer where the particles are stored
     * @param x Array representation three dimensional position of the particle
     * @param n Array representation of quantity of particles in each dimension
     * @param v Array representation of the velocities of the particles in all three dimensions
     * @param h Double representation of distance between the particles
     * @param m Double representation of mass of one particle
     * @param meanVelocity Double representation of mean of the velocity
     * @param sigma sigma value for Lennard Jones
     * @param epsilon epsilon value for Lennard Jones
     * @param type type of the particle
     * @param fz force in z up
     */
    void
    generateMembraneNoBrownian(std::shared_ptr<T> &particles, std::array<double, 3> x, std::array<int, 3> n,
                             std::array<double, 3> v, double h, double m, double sigma = 1, double epsilon = 5, int type = 3, double fz = 1) {
        std::array<double, 3> newX{};
        //iterate over cuboid in each dimension
        std::array<int, 2> index{};
        for (int x_cord = 0; x_cord < n[0]; x_cord++) {
            newX[0] = x[0] + (x_cord * h);
            for (int y_cord = 0; y_cord < n[1]; y_cord++) {
                newX[1] = x[1] + (y_cord * h);
                for (int z_cord = 0; z_cord < n[2]; z_cord++) {
                    newX[2] = x[2] + (z_cord * h);
                    //stores particle in ParticleContainer
                    Particle p = Particle(newX, v, m, sigma, epsilon, type);
                    index[0] = x_cord;
                    index[1] = y_cord;
                    p.setIndex(index);
                    std::array<double, 3> force = {0,0,fz};
                    if((index[0] == 17 && index[1] == 24) || (index[0] == 17 && index[1] == 25 )|| (index[0] == 18 && index[1] == 24 )||( index[0] == 18 && index[1] == 25)){
                        p.setF(force);
                    }
                    particles->addParticle(p);
                }
            }
        }
    }

    /**
     * @brief function to generate sphere from the given parameters
     *
     * @param particles Container where the particles are stored
     * @param center center of the sphere
     * @param r number of particles in each direction
     * @param m mass of each particle
     * @param v velocity of the sphere
     * @param h distance between the particles
     * @param sigma sigma value for Lennard Jones
     * @param epsilon epsilon value for Lennard Jones
     * @param type type of the particle
     */
    void generateSphere(std::shared_ptr<T> &particles, std::array<double, 3> center, int r, double m,
                        std::array<double, 3> v, double h, double sigma = 1, double epsilon = 5, int type=1) {
        double radius = r * h;
        std::array<double, 3> newX{};
        for (int x = -r + 1; x < r; ++x) {
            newX[0] = (x * h);
            for (int y = -r + 1; y < r; ++y) {
                newX[1] = (y * h);
                    double rad_test = ArrayUtils::L2Norm(newX);
                    if (rad_test < radius) {
                        std::array<double, 3> newV = calculateV(v);
                        std::array<double, 3> pos = newX + center;
                        particles->addParticle(Particle(pos, newV, m, sigma, epsilon, type));
                    }

            }
        }

    }

    /**
     * @brief function to generate sphere from the given parameters and the MaxwellBoltzmannDistribution
     * @param particles container where the particles are stored
     * @param center center of the sphere
     * @param v velocity of the sphere
     * @param r number of particles in each direction
     * @param m mass of each particle
     * @param h distance between the particles
     * @param meanVelocity double representation of mean of the velocity
     * @param sigma sigma value for Lennard Jones
     * @param epsilon epsilon value for Lennard Jones
     * @param type type of the particle
     */
    void
    generateSphereBrownian(std::shared_ptr<T> &particles, std::array<double, 3> center, std::array<double, 3> v, int r,
                           double m, double h, double meanVelocity, double sigma = 1, double epsilon = 5, int type=1) {
        double radius = r * h;
        std::array<double, 3> newX{};
        for (int x = -r + 1; x < r; ++x) {
            newX[0] = (x * h);
            for (int y = -r + 1; y < r; ++y) {
                newX[1] = (y * h);

                    double rad_test = ArrayUtils::L2Norm(newX);
                    if (rad_test < radius) {
                        std::array<double, 3> mbV = calculateVWithMean(meanVelocity);
                        std::array<double, 3> newV{};
                        for (int i = 0; i < 3; i++) {
                            newV[i] = v[i] + mbV[i];
                        }
                        std::array<double, 3> pos = newX + center;
                        particles->addParticle(Particle(pos, newV, m, sigma, epsilon, type));
                    }

            }
        }

    }

    /**
     * @brief function to generate sphere from the given parameters
     * @param particles container where the particles are stored
     * @param center center of the sphere
     * @param v velocity of the sphere
     * @param r number of particles in each direction
     * @param m mass of each particle
     * @param h distance between the particles
     * @param sigma sigma value for Lennard Jones
     * @param epsilon epsilon value for Lennard Jones
     * @param type type of the particle
     */
    void generateSphereNoBrownian(std::shared_ptr<T> &particles, std::array<double, 3> center, std::array<double, 3> v, int r,
                                  double m, double h, double sigma = 1, double epsilon = 5, int type=1) {
        double radius = r * h;
        std::array<double, 3> newX{};
        for (int x = -r + 1; x < r; ++x) {
            newX[0] = (x * h);
            for (int y = -r + 1; y < r; ++y) {
                newX[1] = (y * h);

                    double rad_test = ArrayUtils::L2Norm(newX);
                    if (rad_test < radius) {
                        std::array<double, 3> pos = newX + center;
                        particles->addParticle(Particle(pos, v, m, sigma, epsilon, type));
                    }

            }
        }

    }
};

