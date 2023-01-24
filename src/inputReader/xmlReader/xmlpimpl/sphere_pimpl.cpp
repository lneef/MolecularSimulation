//
// Created by dominik on 14.12.22.
//

#include "sphere_pimpl.h"
#include "inputReader/ParticleGenerator.h"
#include "inputReader/xmlReader/LinkedCellStrategy.h"
#include <iostream>

namespace XMLReader {
    void sphere_pimpl::init(std::shared_ptr<LinkedCellStrategy> &lc, std::shared_ptr<Simulation> &simulation) {
        cells = lc;
        sim = simulation;
    }

    void sphere_pimpl::center_x(double x_arg) {
        centre.push(x_arg);
    }

    void sphere_pimpl::center_y(double y_arg) {
        centre.push(y_arg);
    }

    void sphere_pimpl::center_z(double z_arg) {
        centre.push(z_arg);
    }

    void sphere_pimpl::radius(int r) {
        radius_r = r;
    }

    void sphere_pimpl::mesh_width(double h) {
        width = h;
    }

    void sphere_pimpl::mass(double m) {
        mass_r = m;
    }

    void sphere_pimpl::velocity_x(double v_arg) {
        vel.push(v_arg);
    }

    void sphere_pimpl::velocity_y(double v_arg) {
        vel.push(v_arg);
    }

    void sphere_pimpl::velocity_z(double v_arg) {
        vel.push(v_arg);
    }

    void sphere_pimpl::brownianMotion(bool bm) {
        browMot = bm;
    }

    void sphere_pimpl::post_sphere() {
        std::array<double, 3> c{};
        std::array<double, 3> v{};
        size_t len= centre.size();
        for (size_t i = 0; i < len; ++i) {
            c[i] = centre.front();
            v[i] = vel.front();

            centre.pop();
            vel.pop();
        }
        ParticleGenerator<LinkedCellDataStructure> cub{};
        if (!browMot) {
            cub.generateSphereNoBrownian(cells->get(), c, v, radius_r, mass_r, width);
        } else {
            double meanVelocity;
            if (sim->getThermostat() != nullptr) {
                meanVelocity = sqrt(sim->getThermostat()->getTemp() / mass_r);
            } else {
                meanVelocity = 0.1;
            }
            cub.generateSphereBrownian(cells->get(), c, v, radius_r, mass_r, width, meanVelocity, sigma_s, epsilon_s, type_s);
        }

        browMot = true;
        type_s = 1;
        sigma_s = 1;
        epsilon_s = 5;
    }

    void sphere_pimpl::type(int type_arg) {
        type_s = type_arg;
    }

    void sphere_pimpl::sigma(double sigma_arg) {
        sigma_s = sigma_arg;
    }

    void sphere_pimpl::epsilon(double epsilon_arg) {
        epsilon_s = epsilon_arg;
    }
}
