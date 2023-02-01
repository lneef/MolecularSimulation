//
// Created by lukas on 01.12.22.
//

#include "simulation_pimpl.h"
#include "MolSimLogger.h"

namespace XMLReader {
    void simulation_pimpl::t_end(double end_t) {
        sim->setEndTime(end_t);
    }

    void simulation_pimpl::delta_t(double delta_t) {
        sim->setDeltaT(delta_t);
    }

    void simulation_pimpl::domain_size_x(double x) {
        domain.push(x);
    }

    void simulation_pimpl::domain_size_y(double y) {
        domain.push(y);
    }

    void simulation_pimpl::domain_size_z(double z) {
        domain.push(z);
    }

    void simulation_pimpl::cutOff_radius(double cutOff) {
        rCutOff = cutOff;
    }

    void simulation_pimpl::output_name(const std::string &name) {
        sim->setOut_name(name);
    }

    void simulation_pimpl::parallel_mode(const ::std::string & par){
        if(par == "generic"){
            mode = 1;
        }else if(par == "tasking"){
            mode = 2;
        }
    }

    void simulation_pimpl::output_frequency(int f) {
        sim->setOut_frequency(f);
    }

    void simulation_pimpl::l_radius(double rl){
        sim->setForce(std::make_unique<SLennardJones>(rCutOff, rl));
    }

    void simulation_pimpl::post_simulation() {
        std::array<double, 3> dom{};
        size_t i = 0;
        while (!domain.empty()) {
            dom[i++] = domain.front();
            domain.pop();
        }
        auto& chosen = cells -> chose(dim, mode);
        chosen->setSize(rCutOff, dom);
        Simulation::setDim(dim);
        sim->setParticle(chosen);
        MolSimLogger::logInfo("XMLReader: domain=({}, {}, {}), cutoff_radius = {}", dom[0], dom[1], dom[2], rCutOff);
    }


    void simulation_pimpl::init(std::shared_ptr<LinkedCellStrategy> &cells_arg, std::shared_ptr<Simulation> &simu) {
        cells = cells_arg;
        sim = simu;
    }

    void simulation_pimpl::g_gravitation(double g) {
        sim->setForce(std::make_unique<LJGravitation>(g));
        sim->setG(g);
    }

    void simulation_pimpl::dimension(int dim_arg) {
        dim = dim_arg;
    }


}
