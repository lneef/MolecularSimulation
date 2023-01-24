//
// Created by lukas on 01.12.22.
//

#include <iostream>
#include "membrane_pimpl.h"
#include "inputReader/ParticleGenerator.h"
#include "MolSimLogger.h"
#include "inputReader/xmlReader/LinkedCellStrategy.h"

namespace XMLReader {

    void membrane_pimpl::lower_left_x(double x_arg) {
        pos.push(x_arg);
    }

    void membrane_pimpl::lower_left_y(double y_arg) {
        pos.push(y_arg);
    }

    void membrane_pimpl::lower_left_z(double z_arg) {
        pos.push(z_arg);
    }

    void membrane_pimpl::number_x(int n_xarg) {
        num.push(n_xarg);
    }

    void membrane_pimpl::number_y(int n_yarg) {
        num.push(n_yarg);
    }

    void membrane_pimpl::number_z(int n_zarg) {
        num.push(n_zarg);
    }

    void membrane_pimpl::mesh_width(double dist) {
        width = dist;
    }

    void membrane_pimpl::mass(double m_arg) {
        m = m_arg;
    }

    void membrane_pimpl::velocity_x(double v_arg) {
        vel.push(v_arg);
    }

    void membrane_pimpl::velocity_y(double v_arg) {
        vel.push(v_arg);
    }

    void membrane_pimpl::velocity_z(double v_arg) {
        vel.push(v_arg);
    }

    void membrane_pimpl::brownianMotion(bool bm) {
        browMot = bm;
    }

    void membrane_pimpl::fz_up(double f_arg) {
        f = f_arg;
    }

    void membrane_pimpl::stiffness_const(double stiff_const_arg) {
        stiff_const = stiff_const_arg;
    }

    void membrane_pimpl::bond_length(double bond_len_arg) {
        bond_len = bond_len_arg;
    }

    void membrane_pimpl::post_membrane() {
        std::array<double, 3> x{};
        std::array<double, 3> v{};
        std::array<int, 3> n{ 1, 1, 1 };
        size_t i = 0;
        while (!pos.empty()) {
            x[i] = pos.front();
            v[i] = vel.front();
            n[i] = num.front();

            pos.pop();
            vel.pop();
            num.pop();

            ++i;
        }
        ParticleGenerator<LinkedCellDataStructure> cub{};
        if (!browMot) {
            cub.generateMembraneNoBrownian(cells->get(), x, n, v, width, m, sigma_p, epsilon_p, 3, f);
        }
        else {
            double meanVelocity;
            if (sim->getThermostat() != nullptr) {
                meanVelocity = sqrt(sim->getThermostat()->getTemp() / m);
            }
            else {
                meanVelocity = 0.1;
            }
            cub.generateMembraneBrownian(cells->get(), x, n, v, width, m, meanVelocity, sigma_p, epsilon_p, 3, f);
        }

        sim->setIsMembrane(true);
        sim->setF_up(f);

        browMot = true;
        sigma_p = 1;
        epsilon_p = 5;
        f = 1;


        sim->setForce(std::make_unique<MembraneForce>(bond_len, stiff_const));
    }


    void membrane_pimpl::sigma(double sigma_arg) {
        sigma_p = sigma_arg;
    }

    void membrane_pimpl::epsilon(double epsilon_arg) {
        epsilon_p = epsilon_arg;
    }

    void membrane_pimpl::init(std::shared_ptr<LinkedCellStrategy> &lc, std::shared_ptr<Simulation> &sim_arg) {
        cells = lc;
        sim = sim_arg;
    }

}