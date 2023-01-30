//
// Created by lukas on 18.12.22.
//

#include <numeric>
#include <complex>
#include <iostream>
#include "ArrayUtils.h"
#include "Thermostat.h"


void Thermostat::applyThermostat(std::shared_ptr<Container> &particles) {

    double temperature = measureTemp(particles);
    if(temperature == temp_target)
        return;

    double temp_diff = temp_target - temperature;
    if (std::abs(temp_diff) > temp_delta)
        temp_diff = temp_diff > 0 ? temp_delta : -temp_delta;
    double temperature_new = temperature + temp_diff;
    double scaling = std::sqrt(temperature_new / temperature);

    particles->apply([&scaling](Particle &p) {
        std::array<double, 3> newV = scaling * p.getV();
        p.setV(newV);
    });

}

double Thermostat::scalarProduct(Particle &p) {
    auto &v = p.getV();
    return std::accumulate(v.cbegin(), v.cend(), 0., [](double acc, double entry) {
        return acc + entry * entry;
    });
}

void Thermostat::setTarget(double target) {
    temp_target = target;
}

void Thermostat::setDelta(double delta) {
    temp_delta = delta;
}

Thermostat::Thermostat(double current) : temp{current}{}

double Thermostat::getTemp() const{
    return temp;
}

double Thermostat::measureTemp(std::shared_ptr<Container> &particles) {
    double energy = 0.;
    particles->applyPar([&energy](Particle &p) {

#pragma omp atomic
        energy = energy + p.getM() * scalarProduct(p) / 2;
    });
    double temperature = (energy * 2) / (2 * particles->size());
    temp = temperature;

    return temperature;
};
