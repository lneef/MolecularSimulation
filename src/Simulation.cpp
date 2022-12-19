//
// Created by lukas on 08.11.22.
//
#include "Simulation.h"
#include "MolSimLogger.h"
#include "utils/ArrayUtils.h"
#include "outputWriter/XYZWriter.h"
#include <iostream>
#include <iomanip>
void Simulation::calculateX() {
    particles->applyX([this](Particle& p) {
        const std::array<double, 3>& tempV{ p.getV() };
    const std::array<double, 3>& tempF{ p.getF() };
    const std::array<double, 3>& tempX{ p.getX() };

    std::array<double, 3> newX{};

    for (int i = 0; i < 3; i++)
        //Velocity-Störmer-Verlet-Algorithm
        newX[i] = tempX[i] + delta_t * tempV[i] + pow(delta_t, 2) * tempF[i] / (2 * p.getM());

    p.setX(newX);
        });
}

void Simulation::calculateV() {
    particles->apply([this](Particle& p) {
        const std::array<double, 3>& tempV{ p.getV() };
    const std::array<double, 3>& tempOldF{ p.getOldF() };
    const std::array<double, 3>& tempF{ p.getF() };

    std::array<double, 3> newV{};
    for (int i = 0; i < 3; i++) {
        //Velocity-Störmer-Verlet-Algorithm
        newV[i] = tempV[i] + delta_t * (tempOldF[i] + tempF[i]) / (2 * p.getM());
    }

    p.setV(newV);
        });
}


void Simulation::run() {
    //flush logger in benchmark mode
#ifdef BENCHMARK
    MolSimLogger::logger()->flush();
#endif

    auto start = std::chrono::high_resolution_clock::now();

    double current_time = start_time;

    int iteration = 0;

    force->calculateF(particles);

    while (current_time < end_time) {

        calculateX();

        SPDLOG_LOGGER_INFO(MolSimLogger::logger(), "Position of particles calculated for iteration {} ", iteration);
  
        force->calculateF(particles);

        SPDLOG_LOGGER_INFO(MolSimLogger::logger(), "Force on particles calculated for iteration {}", iteration);


        calculateV();
        SPDLOG_LOGGER_INFO(MolSimLogger::logger(), "Velocities of particles calculated for iteration {}", iteration);

        iteration++;

        if (thermostat != nullptr) {
            if (iteration % n_thermostat == 0) {
                thermostat->applyThermostat(particles);
            }
            SPDLOG_LOGGER_INFO(MolSimLogger::logger(), "Temperature: {}", thermostat->measureTemp(particles));
        }


#ifndef BENCHMARK
        if (iteration % out_frequency == 0) {
            writer->plotParticles(particles, out_name, iteration);
        }

        MolSimLogger::logInfo("Itertation {} finished. ", iteration);
#endif
        current_time += delta_t;
    }


    auto stop = std::chrono::high_resolution_clock::now();
    auto difference = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
    std::cout << difference.count() << "ms" << std::endl;

}

Simulation::Simulation(std::shared_ptr<Container>& particles, double delta_t, double end_time,
    std::unique_ptr<outputWriter::FileWriter>& writer, std::unique_ptr<Force>& force) {
    this->delta_t = delta_t;
    this->end_time = end_time;
    this->writer = std::move(writer);
    this->force = std::move(force);
    this->particles = std::move(particles);

}

void Simulation::setDeltaT(double delta_t_arg) {
    delta_t = delta_t_arg;
}

void Simulation::setEndTime(double end_time_arg) {
    end_time = end_time_arg;
}

void Simulation::setParticle(std::shared_ptr<ParticleContainer>& particles_arg) {
    particles = particles_arg;
}

void Simulation::setParticle(std::shared_ptr<LinkedCellContainer>& particles_arg) {
    particles = particles_arg;
}

Simulation::Simulation(double delta_t_arg, double end_time_arg) {
    delta_t = delta_t_arg;
    end_time = end_time_arg;
}

void Simulation::setForce(std::unique_ptr<Force>& force_arg) {
    force = std::move(force_arg);
}

void Simulation::setOut_frequency(int out_frequency_arg) {
    out_frequency = out_frequency_arg;
}

void Simulation::setWriter(std::unique_ptr<outputWriter::FileWriter>& writer_arg) {
    writer = std::move(writer_arg);
}

void Simulation::setOut_name(const std::string& out_name_arg) {
    out_name = out_name_arg;
}

void Simulation::setN_thermostat(int n_thermostat_arg) {
    n_thermostat = n_thermostat_arg;
}


void Simulation::setThermostat(std::shared_ptr<Thermostat>& thermostat_arg) {
    thermostat = thermostat_arg;
}

const std::shared_ptr<Thermostat>& Simulation::getThermostat() const { return thermostat; }

void Simulation::setForce(std::unique_ptr<LJGravitation>&& force_arg) {
    force = std::move(force_arg);
}

const std::unique_ptr<Force>& Simulation::getForce() const {
    return force;
}

void Simulation::checkpoint(const std::string& filename) {
    std::ofstream file;
    std::stringstream strstr;
    strstr << filename << ".txt";

    file.open(strstr.str().c_str());
    file << " #position                           vilocity                             force                               old_force                         masse     type   sigma     epsilon "
        << std::endl;
    file << particles->size() << std::endl;
    particles->apply([&file](Particle& p) {

        //print the position of each particle
        std::array<double, 3> x = p.getX();
    file.setf(std::ios_base::showpoint);
    for (auto& xi : x) {
        file << std::setw(10) << xi << " ";
    } file << "   ";

    //print the velocity of each particle
    std::array<double, 3> v = p.getV();
    file.setf(std::ios_base::showpoint);
    for (auto& vi : v) {
        file << std::setw(10) << vi << " ";
    }file << "   ";

    //print the force of each particle
    std::array<double, 3> f = p.getF();
    file.setf(std::ios_base::showpoint);
    for (auto& fi : f) {
        file << std::setw(10) << fi << " ";
    }file << "   ";

    //print the old_force of each particle
    std::array<double, 3> old_f = p.getOldF();
    file.setf(std::ios_base::showpoint);
    for (auto& fi : old_f) {
        file << std::setw(10) << fi << " ";
    }file << "   ";

    //print the masse of each particle
    file << p.getM() << "    ";

    //print the type of each particle
    file << p.getType() << "    ";

    //print the sigma of each particle
    file << p.getSigma() << "    ";

    //print the epsilon of each particle
    file << p.getEpsilon() << " ";

    file << std::endl;
        });

    file.close();
};