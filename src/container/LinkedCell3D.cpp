#include <iostream>
#include "LinkedCell3D.h"
#include "MolSimLogger.h"

LinkedCell3D::~LinkedCell3D() = default;

size_t LinkedCell3D::size() {
    size_t len = 0;

    for (size_t i = 1; i < layers.size() - 1; ++i) {
        len += layers[i].size();
    }

    return len;
}

void LinkedCell3D::apply(std::function<void(Particle &)> fun) {
    for (size_t i = 1; i < layers.size() - 1; ++i) {
        layers[i].apply(fun);
    }

}

void LinkedCell3D::addParticle(Particle &&p) {
    auto ind = index(p.getX());
    if (ind >= mesh[2]) {
        return;
    }
    layers[ind].addParticle(p);
}

size_t LinkedCell3D::index(const std::array<double, 3> &pos) noexcept {

    if (pos[2] >= domain[2]) {
        return mesh[2] - 1;
    }

    size_t ind = std::floor(pos[2] / cutoff[2]) + 1;

    return ind;
}

void LinkedCell3D::applyF(std::function<void(Particle &, Particle &)> fun) {
#pragma omp parallel shared(layers)
    {
        preparePeriodic();



#pragma omp for schedule(dynamic, 1) nowait
            for (size_t i = 0; i < layerSize; ++i) {
                layers[0][i].apply([this, &fun, i](Particle &p) {
                    forceThreeD(p, i, 0, fun);
                });

            }

#pragma omp for schedule(guided, mesh[0] + 1) collapse(2) nowait
        for (std::size_t j = 1; j < layers.size() - 1; ++j) {
            for (size_t i = 0; i < layerSize; ++i) {
                layers[j].forceTwoD(layers[j][i], i, fun);
                layers[j][i].apply([this, &fun, j, i](Particle &p) {
                    forceThreeD(p, i, j, fun);
                });
            }

        }

#pragma omp for schedule(dynamic, 1)
        for(size_t i = 1; i< layers.size() - 1; ++i){
            layers[i].applyFBoundary(fun);
        }
    }

    applyFBoundary(fun);
}

void
LinkedCell3D::forceThreeD(Particle &p, size_t ind2D, size_t ind3D,
                          std::function<void(Particle &, Particle &)> fun) {

    auto partial = [&p, &fun](Particle &p1) {
        fun(p, p1);
    };
    auto &next = layers[ind3D + 1];
    int width = static_cast<int>(mesh[0]);

    //calculate interaction with all three cells in the next layer below, at the same y-coodinate and above the current cell
    for (int j = -width; j <= width; j += width) {

        if (ind2D + j < 0 || ind2D + j >= next.cells.size())
            continue;


        next[ind2D + j].apply(partial);


        if (ind2D % mesh[0] == mesh[0] - 1) {
            if (ind2D + j - 1 < 0)
                continue;
            next[ind2D + j - 1].apply(partial);
        } else if (ind2D % mesh[0] == 0) {
            if (ind2D + j + 1 > next.cells.size())
                continue;
            next[ind2D + j + 1].apply(partial);
        } else {
            next[ind2D + j + 1].apply(partial);
            next[ind2D + j - 1].apply(partial);

        }

    }
}

void LinkedCell3D::applyX(std::function<void(Particle &)> fun) {
    applyPar(fun);

    update();

    clearHalo();
}

void LinkedCell3D::update() {
       //vectors to store particles collected by the worker threads
    std::vector<std::vector<Particle>> temp;
    std::vector<std::vector<size_t>> temp_ind;
    std::vector<std::vector<size_t>> temp_ind3;
    temp.resize(layerSize);
    temp_ind3.resize(layerSize);
    temp_ind.resize(layerSize);
    //worker threads collect particles that need to be relocated
    for (size_t i = 1; i < layers.size() - 1; ++i) {

        #pragma omp parallel for schedule(dynamic, 1)
        for (size_t j = mesh[0] + 1; j < layerSize - mesh[0] - 1; ++j) {
            if(j % mesh[0] == 0 || j % mesh[0] == mesh[0] - 1)
                continue;
            for (auto it = layers[i][j].begin(); it != layers[i][j].end();) {
                auto &p = *it;
                size_t ind = layers[i].index(p);
                size_t ind3D = index(p.getX());


                if (ind == j && ind3D == i) {
                    ++it;
                    continue;
                }

                temp[j].push_back(*it);
                temp_ind[j].push_back(ind);
                temp_ind3[j].push_back(ind3D);

                it = layers[i][j].remove(it);


            }
        }

        //relocation is done by the master thread
        for(size_t l = 0; l < layerSize; ++l){
            for(size_t k = 0; k<temp[l].size(); ++k){
                auto &pos = temp[l][k].getX();

                if (pos[0] < -cutoff[0] || pos[0] >= domain[0] + cutoff[0] || pos[1] < -cutoff[1] ||
                    pos[1] > domain[1] + cutoff[1]) {
                    continue;
                }


                LinkedCell3D::update(temp[l][k], temp_ind3[l][k], temp_ind[l][k]);
            }
            temp[l].clear();
            temp_ind[l].clear();
            temp_ind3[l].clear();
        }


    }

}


void LinkedCell3D::clearHalo() {
#pragma omp parallel for schedule(static, 1)
    for (size_t i = 1; i < layers.size() - 1; ++i) {
        layers[i].clearHalo();
    }

    for (auto &cell: layers[0].cells)
        cell.clear();


    for (auto &cell: layers[layers.size() - 1].cells)
        cell.clear();

}

void LinkedCell3D::preparePeriodic() {

#pragma omp for schedule(dynamic, 1)
    for (size_t i = 1; i < layers.size() - 1; ++i) {
        layers[i].updatePeriodic();
    }

#pragma omp single
    {
        //mirror inner cells
        if (containsPeriodic(Boundary::BACK) || containsPeriodic(Boundary::FRONT)) {
            frontBackBoundary(-domain[2], layers.size() - 2, 0);
            frontBackBoundary(domain[2], 1, layers.size() - 1);
        }

        //mirror left and right boundary
        if (containsPeriodic(Boundary::BACK) || containsPeriodic(Boundary::FRONT) ||
            containsPeriodic(Boundary::LEFT) ||
            containsPeriodic(Boundary::RIGHT)) {
            mirrorHorizontal(domain[2], 1, layers[layers.size() - 1]);
            mirrorHorizontal(-domain[2], layers.size() - 2, layers[0]);
        }

        //mirror top and bottom boundary
        if (containsPeriodic(Boundary::BACK) || containsPeriodic(Boundary::FRONT) ||
            containsPeriodic(Boundary::BOTTOM) ||
            containsPeriodic(Boundary::TOP)) {
            mirrorVertical(domain[2], 1, layers[layers.size() - 1]);
            mirrorVertical(-domain[2], layers.size() - 2, layers[0]);
        }

        //mirror corner cells
        if (!periodic.empty()) {
            mirrorDiagonal(domain[2], 1, layers[layers.size() - 1]);
            mirrorDiagonal(-domain[2], layers.size() - 2, layers[0]);

        }
    }

}

void LinkedCell3D::frontBackBoundary(double to_add, size_t ind, size_t oth) {

    auto &lc = layers[ind];
    auto &counter = layers[oth];

    size_t i = LinkedCellContainer::mesh[0] + 1;
    for (; i < lc.cells.size() - mesh[0] - 1; ++i) {

        if (i % mesh[0] == 0 || i % mesh[0] == mesh[0] - 1)
            continue;

        for (auto &p: lc[i]) {
            auto newP = p.getX();
            newP[2] += to_add;
            counter.simpleAdd(Particle(newP, p.getV(), p.getM(), p.getSigma(), p.getEpsilon(), p.getType(), true));
        }
    }


}

bool LinkedCell3D::side(size_t ind3D) {
    return (ind3D == 0 && containsPeriodic(Boundary::FRONT)) ||
           (ind3D == layers.size() - 1 && containsPeriodic(Boundary::BACK));

}

void LinkedCell3D::update(Particle &particle, size_t ind3D, size_t ind) {

    //move particle if periodic boundaries for fornt and back are specified
    if (side(ind3D)) {
        ind3D = updatePeriodic(particle, ind3D);

        //move particle if it front and back is not specified, but it is contained in a boundary cell for which a 2D periodic boundary is specified
    }else if((ind3D == 0 || ind3D == layers.size() - 1) && layers[ind3D].side(ind)){
        ind3D = updatePeriodic(particle, ind3D);
    }

    if (layers[ind3D].side(ind)) {
        ind = layers[ind3D].mirror(particle, ind);

    }
    if (ind >= layerSize)
        return;

    layers[ind3D][ind].addParticle(particle);

}

size_t LinkedCell3D::updatePeriodic(Particle &p, size_t ind3D) {
    auto pos = p.getX();
    if (ind3D == 0) {
        pos[2] += domain[2];
        ind3D = layers.size() - 2;
    } else {
        pos[2] -= domain[2];
        ind3D = 1;
    }
    p.setX(pos);
    return ind3D;
}

void LinkedCell3D::setSize(double cutOff_arg, std::array<double, 3> &domain_arg) {
    domain = domain_arg;
    for (size_t i = 0; i < 3; ++i) {
        mesh[i] = std::floor(std::abs(domain_arg[i]) / cutOff_arg);
        cutoff[i] = domain[i] / static_cast<double>(mesh[i]);
        mesh[i] +=2;
    }

    layers.resize(mesh[2]);
    for (auto &layer: layers) {
        layer.setUp();
    }

    layerSize = mesh[0] * mesh[1];

}

void LinkedCell3D::addParticle(Particle &p) {
    auto &pos = p.getX();
    size_t ind3D = index(pos);
    layers[ind3D].addParticle(p);

}

std::array<double, 3> &LinkedCell3D::getDomain() {
    return domain;
}

void LinkedCell3D::setDomain(std::array<double, 3> &domain_arg) {
    domain = domain_arg;
}

LinkedCellContainer &LinkedCell3D::operator[](size_t i) {
    return layers[i];
}

void LinkedCell3D::mirrorHorizontal(double dis, size_t i, LinkedCellContainer &counter) {

    for (auto &cell: layers[i].left_boundary) {
        std::array<double, 3> newP{domain[0], 0, dis};
        for (auto &p: cell.get()) {
            counter.simpleAdd(
                    Particle(p.getX() + newP, p.getV(), p.getM(), p.getSigma(), p.getEpsilon(), p.getType(), true));
        }
    }

    for (auto &cell: layers[i].right_boundary) {
        std::array<double, 3> newP{-domain[0], 0, dis};
        for (auto &p: cell.get()) {
            counter.simpleAdd(
                    Particle(p.getX() + newP, p.getV(), p.getM(), p.getSigma(), p.getEpsilon(), p.getType(), true));
        }
    }


}


void LinkedCell3D::mirrorVertical(double dis, size_t i, LinkedCellContainer &counter) {

    for (auto &cell: layers[i].bottom_boundary) {
        std::array<double, 3> newP{0, domain[1], dis};
        for (auto &p: cell.get()) {
            counter.simpleAdd(
                    Particle(p.getX() + newP, p.getV(), p.getM(), p.getSigma(), p.getEpsilon(), p.getType(), true));
        }
    }

    for (auto &cell: layers[i].top_boundary) {
        std::array<double, 3> newP{0, -domain[1], dis};
        for (auto &p: cell.get()) {
            counter.simpleAdd(
                    Particle(p.getX() + newP, p.getV(), p.getM(), p.getSigma(), p.getEpsilon(), p.getType(), true));
        }
    }
}

void LinkedCell3D::mirrorDiagonal(double dis, size_t i, LinkedCellContainer &counter) {
    //mirror particles from corner cells to other corner cells

    for (auto &p: layers[i].bottom_boundary[0].get()) {
        std::array<double, 3> to_add{domain[0], domain[1], dis};
        counter.simpleAdd(
                Particle(p.getX() + to_add, p.getV(), p.getM(), p.getSigma(), p.getEpsilon(), p.getType(), true));
    }

    for (auto &p: layers[i].top_boundary[0].get()) {
        std::array<double, 3> to_add{domain[0], -domain[1], dis};
        counter.simpleAdd(
                Particle(p.getX() + to_add, p.getV(), p.getM(), p.getSigma(), p.getEpsilon(), p.getType(), true));
    }

    for (auto &p: layers[i].right_boundary[0].get()) {
        std::array<double, 3> to_add{-domain[0], domain[1], dis};
        counter.simpleAdd(
                Particle(p.getX() + to_add, p.getV(), p.getM(), p.getSigma(), p.getEpsilon(), p.getType(), true));
    }

    size_t len = layers[i].top_boundary.size();
    for (auto &p: layers[i].top_boundary[len - 1].get()) {
        std::array<double, 3> to_add{-domain[0], -domain[1], dis};
        counter.simpleAdd(
                Particle(p.getX() + to_add, p.getV(), p.getM(), p.getSigma(), p.getEpsilon(), p.getType(), true));
    }

}


void LinkedCell3D::applyFBoundary(std::function<void(Particle &, Particle &)> fun) {

    if (conditions.find(Boundary::FRONT) != conditions.end()) {
        auto &cond = conditions.at(Boundary::FRONT);
        for (auto &cell: layers[1].cells) {
            cell.apply([&fun, &cond](Particle &p) {
                if (cond.check(p)) {
                    cond.apply(p, fun);
                }
            });
        }
    }

    if (conditions.find(Boundary::BACK) != conditions.end()) {
        auto &cond = conditions.at(Boundary::BACK);
        for (auto &cell: layers[mesh[0] - 1].cells) {
            cell.apply([&fun, &cond](Particle &p) {
                if (cond.check(p)) {
                    cond.apply(p, fun);
                }
            });
        }
    }
}

void LinkedCell3D::applyPar(std::function<void(Particle &)> fun) {
#pragma omp parallel for schedule(dynamic, 1) collapse(2)
    for (size_t i = 1; i < layers.size() - 1; ++i) {
        for (size_t j = mesh[0] + 1; j < layerSize - mesh[0] - 1; ++j){
            if(j% mesh[0] == 0 || j % mesh[0] == mesh[0] - 1)
                continue;

            layers[i][j].apply(fun);
        }
    }
}

LinkedCell3D::LinkedCell3D() =
default;
