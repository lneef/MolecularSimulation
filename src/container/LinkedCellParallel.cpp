#include "LinkedCellParallel.h"

void LinkedCellParallel::applyF(std::function<void(Particle &, Particle &)> fun) {
#pragma omp parallel shared(layers)
    {
        preparePeriodic();
#pragma omp for
        for (size_t i = 0; i < layerSize; ++i) {
            layers[0][i].apply([this, &fun, i](Particle &p) {
                forceThreeD(p, i, 0, fun);
            });

        }

#pragma omp single
        {

            for (size_t z = 1; z < mesh[2] - 1; z += 3) {
                for (size_t y = 0; y < mesh[1]; y += 3) {
                    for (size_t x = 0; x < mesh[0]; x += 3) {
                        size_t begin = y * mesh[0] + x;
                        size_t end_x = x + 3 < mesh[0] ? 3 : mesh[0] - x;
                        end_x+=begin;
                        size_t end_z = z + 3 < mesh[2] - 1 ? z + 3 : mesh[2] - 1;
                        size_t end_y = y + 3 < mesh[1] ? 3 : mesh[1] - y;


#pragma omp task shared(layers)
                        {
                            for (size_t z_i = z; z_i < end_z; ++z_i) {
                                for (size_t y_i = 0; y_i < end_y * mesh[0]; y_i += mesh[0]) {
                                    for (size_t x_i = begin; x_i < end_x; ++x_i) {
                                        if(x_i + y_i > layerSize)
                                            break;
                                        layers[z_i].forceTwoD(layers[z_i][x_i + y_i], x_i + y_i, fun);
                                        layers[z_i][x_i + y_i].apply([this, &fun, z_i, x_i, y_i](Particle &p) {
                                            forceThreeD(p, x_i + y_i, z_i, fun);
                                        });
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }


#pragma omp for
        for(size_t i=1; i<layers.size() -1 ; ++i){
            layers[i].applyFBoundary(fun);

        }
    }

    applyFBoundary(fun);
}

void LinkedCellParallel::applyX(std::function<void(Particle &)> fun) {
    applyPar(fun);

    update();

    clearHalo();

}

void LinkedCellParallel::applyPar(std::function<void(Particle &)> fun) {
#pragma omp parallel for collapse(2)
    for(size_t i = 1 ; i<layers.size() - 1; ++i){
        for(size_t j = mesh[0] + 1; j < layerSize - mesh[0] - 1; ++j){
            if(j% mesh[0] == 0 || j % mesh[0] == mesh[0] - 1){
                continue;
            }

            layers[i][j].apply(fun);
        }
    }
}

LinkedCellParallel::~LinkedCellParallel() = default;

LinkedCellParallel::LinkedCellParallel() = default;


size_t LinkedCellParallel::size() {
    size_t len = 0;

    for (size_t i = 1; i < layers.size() - 1; ++i) {
        len += layers[i].size();
    }

    return len;
}

void LinkedCellParallel::apply(std::function<void(Particle &)> fun) {
    for (size_t i = 1; i < layers.size() - 1; ++i) {
        layers[i].apply(fun);
    }

}

void LinkedCellParallel::addParticle(Particle &&p) {
    auto ind = index(p.getX());
    if (ind >= mesh[2]) {
        return;
    }
    layers[ind].addParticle(p);
}



void LinkedCellParallel::setSize(double cutOff_arg, std::array<double, 3> &domain_arg) {
    domain = domain_arg;
    for (size_t i = 0; i < 3; ++i) {
        mesh[i] = std::floor(std::abs(domain_arg[i]) / cutOff_arg);
        cutoff[i] = domain[i] / static_cast<double>(mesh[i]);
        mesh[i] +=2;
    }

    LinkedCellContainer::setDomain(domain);
    LinkedCellContainer::setMesh(mesh);

    layers.resize(mesh[2]);
    for (auto &layer: layers) {
        layer.setUp();
    }

    layerSize = mesh[0] * mesh[1];

}

void LinkedCellParallel::addParticle(Particle &p) {
    auto &pos = p.getX();
    size_t ind3D = index(pos);
    layers[ind3D].addParticle(p);

}

