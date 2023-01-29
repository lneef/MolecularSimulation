//
// Created by lukas on 01.12.22.
//

#include <iostream>
#include "cuboid_pimpl.h"
#include "inputReader/ParticleGenerator.h"
#include "MolSimLogger.h"
#include "inputReader/xmlReader/LinkedCellStrategy.h"

namespace XMLReader {

    void cuboid_pimpl::lower_left_x(double x_arg) {
        pos.push(x_arg);
    }

    void cuboid_pimpl::lower_left_y(double y_arg) {
        pos.push(y_arg);
    }

    void cuboid_pimpl::lower_left_z(double z_arg) {
        pos.push(z_arg);
    }

    void cuboid_pimpl::number_x(int n_xarg) {
        num.push(n_xarg);
    }

    void cuboid_pimpl::number_y(int n_yarg) {
        num.push(n_yarg);
    }

    void cuboid_pimpl::number_z(int n_zarg) {
        num.push(n_zarg);
    }

    void cuboid_pimpl::mesh_width(double dist) {
        width = dist;
    }

    void cuboid_pimpl::mass(double m_arg) {
        m = m_arg;
    }

    void cuboid_pimpl::velocity_x(double v_arg) {
        vel.push(v_arg);
    }

    void cuboid_pimpl::velocity_y(double v_arg) {
        vel.push(v_arg);
    }

    void cuboid_pimpl::velocity_z(double v_arg) {
        vel.push(v_arg);
    }

    void cuboid_pimpl::brownianMotion(bool bm) {
        browMot = bm;
    }

    void cuboid_pimpl::post_cuboid() {
        std::array<double, 3> x{};
        std::array<double, 3> v{};
        std::array<int, 3> n{1, 1, 1};
        size_t i = 0;
        while (!pos.empty()){
            x[i] = pos.front();
            v[i] = vel.front();
            n[i] = num.front();

            pos.pop();
            vel.pop();
            num.pop();

            ++i;
        }
        ParticleGenerator<LinkedCellDataStructure> cub{Simulation::getDim()};
        if (!browMot) {
            cub.generateCuboidNoBrownian(cells->get(), x, n, v, width, m, sigma_p, epsilon_p, type_p);
        } else {
            double meanVelocity;
            if (sim->getThermostat() != nullptr) {
                meanVelocity = sqrt(sim->getThermostat()->getTemp()/m);
            } else {
                meanVelocity = 0.1;
            }
            cub.generateCuboidBrownian(cells->get(), x, n, v, width, m, meanVelocity, sigma_p, epsilon_p, type_p);
        }

        browMot = true;
        type_p = 1;
        sigma_p = 1;
        epsilon_p = 5;

    }

    void cuboid_pimpl::init(std::shared_ptr<LinkedCellStrategy> &lc, std::shared_ptr<Simulation> &sim_arg) {
        cells = lc;
        sim = sim_arg;
    }

    void cuboid_pimpl::type(int type_arg) {
        type_p = type_arg;
    }

    void cuboid_pimpl::sigma(double sigma_arg) {
        sigma_p = sigma_arg;
    }

    void cuboid_pimpl::epsilon(double epsilon_arg) {
        epsilon_p = epsilon_arg;
    }

}