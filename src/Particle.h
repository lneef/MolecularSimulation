/*
 * Particle.h
 *
 *  Created on: 23.02.2010
 *      Author: eckhardw
 */

#pragma once
#include <omp.h>
#include <array>
#include <string>
#include <memory>
class Particle {

private:
    bool ghost;
#ifdef _OPENMP
    omp_lock_t par_lock;
#endif
    /**
     * Position of the particle
     */
    std::array<double, 3> x;

    /**
     * Velocity of the particle
     */
    std::array<double, 3> v;

    /**
     * Force effective on this particle
     */
    std::array<double, 3> f;

    /**
     * Force which was effective on this particle
     */
    std::array<double, 3> old_f;

    /**
     * the x/y-indices of the particle
     */
    std::array<int, 2> membrane_index;

    /**
    * Mass of this particle
    */
    double m;

    /**
     * Type of the particle. Use it for whatever you want (e.g. to separate
     * molecules belonging to different bodies, matters, and so on)
     */
    int type;

    /**
     * @brief sigma value for the Lennard Jones potential
     */
    double sigma;

    /**
     * @brief epsilon value for the Lennard Jones potential
     */
    double epsilon;

    std::array<double, 3> old_x;

    std::array<int, 3> warp{};
    std::array<double, 3> last_x;

public:
    std::array<int, 3>& getWarp();
    std::array<double, 3>& getLast();
    void setWarp(std::array<int, 3>& warp_arg);

    explicit Particle(int type = 0);

    /**
     * @brief copy constructor of the class particle
     * @param other Particle to be copied
     */
    Particle(const Particle& other);

    void subtractFromF(std::array<double, 3>& to_sub);
    void addToF(std::array<double, 3>& to_add);

    /**
     * @brief constructor of the class particle
     * @param x_arg position of the particles
     * @param v_arg velocity of the particles
     * @param m_arg mass of the particle
     * @param sigma_arg sigma value for the Lennard Jones potential (default: 1)
     * @param epsilon_arg epsilon value for the Lennard Jones potential (default: 5)
     * @param type type of the particle (default: 1)
     */
    Particle(
        // for visualization, we need always 3 coordinates
        // -> in case of 2d, we use only the first and the second
        std::array<double, 3> x_arg, std::array<double, 3> v_arg, double m_arg, double sigma_arg = 1, double epsilon_arg = 5,
        int type = 0, bool ghost_arg = false);

    virtual ~Particle();

    [[nodiscard]] const std::array<double, 3>& getX() const;

    [[nodiscard]] const std::array<double, 3>& getV() const;

    [[nodiscard]] const std::array<double, 3>& getF() const;

    [[nodiscard]] const std::array<double, 3>& getOldF() const;

    [[nodiscard]] const std::array<int, 2>& getIndex() const;

    [[nodiscard]] double getM() const;

    [[nodiscard]] int getType() const;

    [[nodiscard]] const std::array<double, 3>& getOldX() const;

    bool operator==(Particle& other);

    [[nodiscard]] std::string toString() const;

    /**
     * @brief getter for field sigma
     * @return sigma value for Lennard Jones potential of the Particle
     */
    [[nodiscard]] double getSigma() const;

    /**
   * @brief getter for field epsilon
   * @return epsilon value for Lennard Jones potential of the Particle
   */
    [[nodiscard]] double getEpsilon() const;

    /**
     * @brief Setter for the effective force acting on the particle. Before f is set to the new value, old_f is set to the old value of f.
     *
     * @param f Array representing three dimensional vector of forces acting on the particle
    */
    void setF(const std::array<double, 3>& f);

    /**
     * @brief Setter for the position of the particle.
     *
     * @param x Array representing three dimensional position of the particle
    */

    void setX(const std::array<double, 3>& x);

    /**
     * @brief Setter for the velocity of the particle
     *
     * @param v three dimensional vector of the velocities of the particle in all three dimensions
    */
    void setV(const std::array<double, 3>& v);

    /**
     * @brief function which sets f to a new value, it does not update old_f
     * @param f Array representing three dimensional vector of the forces acting on the particle
     *
     * @warning old_f is set to old value of f
     */
    void updateF(const std::array<double, 3>& f);

    /**
     * @brief Setter for mass of the particle
     * @param m double representing new mass of the particle
     */
    void setM(double m);

    /**
     * @brief Setter for old_f of the particle
     * @param oldf Array representing three dimensional vector of forces that were effective on the particle
     */
    void setOldF(const std::array<double, 3>& oldf);

    /**
     * @brief Setter for membrane_index of the particle
     * @param index Array representing two dimensional indices of the membrane-particle
     */
    void setIndex(const std::array<int, 2>& index);

    void setOldX(const std::array<double, 3>& oldx_arg);

    /**
     * @brief to detect if the another membrane-particle is the direct neighbor  
     * @param another membrane-particle
     */
    bool ifDirectNeighbor(const Particle& p);

    /**
     * @brief to detect if the another membrane-particle is the diagonal neighbor  
     * @param another membrane-particle
     */
    bool ifDiagonalNeighbor(const Particle& p);

    /**
     * @brief compares two values of type double for equality
     * @param d1 first value of type double
     * @param d2 second value of type double
     * @return true if both values are approximately equal, false otherwise
     */
    static bool comp(double d1, double d2);
};

std::ostream& operator<<(std::ostream& stream, Particle& p);
