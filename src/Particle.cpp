/*
 * Particle.cpp
 *
 *  Created on: 23.02.2010
 *      Author: eckhardw
 */

#include "Particle.h"

#include <iostream>
#include "utils/ArrayUtils.h"
#include "MolSimLogger.h"

Particle::Particle(int type_arg) {
    type = type_arg;
    MolSimLogger::logTrace("Particle generated!");
    f = { 0., 0., 0. };
    old_f = { 0., 0., 0. };
    ghost = false;
#ifdef __OPENMP
    omp_init_lock(par_lock);
#endif
}

Particle::Particle(const Particle& other) {
    x = other.x;
    v = other.v;
    f = other.f;
    old_f = other.old_f;
    m = other.m;
    type = other.type;
    sigma = other.sigma;
    epsilon = other.epsilon;
    membrane_index = other.membrane_index;
    ghost = other.ghost;
#ifdef __OPENMP
    omp_init_lock(par_lock);
#endif
    MolSimLogger::logTrace("Particle generated by copy!");
}

// Todo: maybe use initializater list instead of copy?
Particle::Particle(std::array<double, 3> x_arg, std::array<double, 3> v_arg,
    double m_arg, double sigma_arg, double epsilon_arg, int type_arg, bool ghost_arg) {
    x = x_arg;
    v = v_arg;
    m = m_arg;
    type = type_arg;
    f = { 0., 0., 0. };
    old_f = { 0., 0., 0. };
    sigma = sigma_arg;
    epsilon = epsilon_arg;
    ghost = ghost_arg;
#ifdef __OPENMP
    omp_init_lock(par_lock);
#endif
    MolSimLogger::logTrace("Particle generated!");
}


Particle::~Particle() { MolSimLogger::logTrace("Particle destructed!"); }

const std::array<double, 3>& Particle::getX() const { return x; }

const std::array<double, 3>& Particle::getV() const { return v; }

const std::array<double, 3>& Particle::getF() const { return f; }

const std::array<double, 3>& Particle::getOldF() const { return old_f; }

const std::array<int, 2>& Particle::getIndex() const { return membrane_index; }

void Particle::updateF(const std::array<double, 3>& f) {
    this->old_f = this->f;
    this->f = f;
}

double Particle::getM() const { return m; }

int Particle::getType() const { return type; }

std::string Particle::toString() const {
    std::stringstream stream;
    stream << "Particle: X:" << x << " v: " << v << " f: " << f
        << " old_f: " << old_f << " type: " << type;
    return stream.str();
}

void Particle::setF(const std::array<double, 3>& f) {
    this->f = f;
}

void Particle::setV(const std::array<double, 3>& v) {
    this->v = v;
}

void Particle::setX(const std::array<double, 3>& x) {
    this->x = x;
}

bool Particle::operator==(Particle& other) {
    return (x == other.x) and (v == other.v) and (f == other.f) and
        (type == other.type) and (m == other.m) and (old_f == other.old_f);
}

std::ostream& operator<<(std::ostream& stream, Particle& p) {
    stream << p.toString();
    return stream;
}

void Particle::setM(const double m) {
    this->m = m;
}

void Particle::setOldF(const std::array<double, 3>& oldf) {
    this->old_f = oldf;
}


void Particle::setIndex(const std::array<int, 2>& index) {
    this->membrane_index = index;
}


bool Particle::comp(double d1, double d2) {
    return std::abs(d1 - d2) < std::numeric_limits<double>::epsilon();
}

double Particle::getSigma() const {
    return sigma;
}

double Particle::getEpsilon() const {
    return epsilon;
}

bool Particle::ifDirectNeighbor(const Particle& p) {
    int x_diff = this->membrane_index.at(0) - p.getIndex().at(0);
    int y_diff = this->membrane_index.at(1) - p.getIndex().at(1);
    return ((x_diff == 0 && (y_diff == 1 || y_diff == -1)) || (y_diff == 0 && (x_diff == 1 || x_diff == -1)));
}


bool Particle::ifDiagonalNeighbor(const Particle& p) {
    int x_diff = this->membrane_index.at(0) - p.getIndex().at(0);
    int y_diff = this->membrane_index.at(1) - p.getIndex().at(1);
    return ((x_diff == 1 && (y_diff == 1 || y_diff == -1)) || (x_diff == -1 && (y_diff == 1 || y_diff == -1)));

}

void Particle::subtractFromF(std::array<double, 3> &to_sub) {
#ifdef __OPENMP
    omp_set_lock(par_lock);
#endif
    f[0] -= to_sub[0];
    f[1] -= to_sub[1];
    f[2] -= to_sub[2];
#ifdef __OPENMP
    omp_unset_lock(par_lock);
#endif
}

void Particle::addToF(std::array<double, 3> &to_add) {
#ifdef __OPENMP
    omp_set_lock(par_lock);
#endif
        f[0] += to_add[0];
        f[1] += to_add[1];
        f[2] += to_add[2];
#ifdef __OPENMP
    omp_unset_lock(par_lock);
#endif
}
