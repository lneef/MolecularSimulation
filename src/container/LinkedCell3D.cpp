//
// Created by lukas on 02.01.23.
//

#include <iostream>
#include "LinkedCell3D.h"
#include "MolSimLogger.h"

LinkedCell3D::~LinkedCell3D() = default;

size_t LinkedCell3D::size() {
    size_t len = 0;

    for (auto &layer: layers)
        len += layer.size();

    return len;
}

void LinkedCell3D::apply(std::function<void(Particle &)> fun) {
    for (auto &layer: layers) {
        layer.apply(fun);
    }

}

void LinkedCell3D::addParticle(Particle &&p) {
    auto ind = index(p.getX());
    layers[ind].addParticle(p);
}

size_t LinkedCell3D::index(const std::array<double, 3> &pos) noexcept {

    if (pos[2] >= domain[2]) {
        return mesh[2] - 1;
    }

    size_t ind = std::floor(pos[2] / cutOff) + cutOff;

    return ind;
}

void LinkedCell3D::applyF(std::function<void(Particle &, Particle &)> fun) {
    size_t ind3d = 0;
    preparePeriodic();
    for (std::size_t j = 0; j < layers.size() - 1; ++j) {
        auto &layer = layers[j];
        for (size_t i = 0; i < layer.cells.size(); ++i) {
            layer.forceTwoD(layer[i], i, fun);
            layer[i].apply([this, &fun, ind3d, i](Particle &p) {
                forceThreeD(p, i, ind3d, fun);
            });
        }
        ind3d++;

    }

}

void
LinkedCell3D::forceThreeD(Particle &p, size_t ind2D, size_t ind3D, std::function<void(Particle &, Particle &)> fun) {

    auto partial = [&p, &fun](Particle &p1) {
        fun(p, p1);
    };
    auto &next = layers[ind3D + 1];
    int width = mesh[0];
    for (int j = -width; j <= width; j += width) {

        if (ind2D + j < 0 || ind2D + j >= next.cells.size())
            continue;


        next[ind2D + j].apply(partial);


        if (ind2D % mesh[0] == mesh[0] - 1) {
            if(ind2D + j - 1 < 0)
                continue;
            next[ind2D + j - 1].apply(partial);
        }else if (ind2D % mesh[0] == 0) {
            if(ind2D + j + 1 > next.cells.size())
                continue;
            next[ind2D + j + 1].apply(partial);
        } else {
            next[ind2D + j + 1].apply(partial);
            next[ind2D + j - 1].apply(partial);

        }

    }
}

void LinkedCell3D::applyX(std::function<void(Particle &)> fun) {
    for (auto &layer: layers) {
        layer.apply(fun);
    }
    update();
}

void LinkedCell3D::update() {
    for (size_t i = 1; i < layers.size() - 1; ++i) {
        for (size_t j = mesh[0] + 1; j < layers[i].cells.size() - mesh[0] - 1; ++j) {
            for (auto it = layers[i][j].begin(); it != layers[i][j].end();) {
                auto &p = *it;
                size_t ind = layers[i].index(p);
                size_t ind3D = index(p.getX());
                auto &pos = p.getX();

                //check that not completely outside
                if (pos[0] < -cutOff || pos[0] >= domain[0] + cutOff || pos[1] < -cutOff ||
                    pos[1] > domain[1] + cutOff) {
                    SPDLOG_LOGGER_INFO(MolSimLogger::logger(), "Particle at position ({}, {}, {}) removed", p.getX()[0],
                                       p.getX()[1], p.getX()[2]);
                    it = layers[i][j].remove(it);
                    continue;
                }

                if (ind == j && ind3D == i) {
                    ++it;
                    continue;
                }

                update(p, ind3D, ind);

                it = layers[i][j].remove(it);

            }

            //skip halo
            if (j % mesh[0] == mesh[0] - 2)
                j += 2;
        }
    }


}


void LinkedCell3D::clearHalo() {
    for (auto &cont: layers)
        cont.clearHalo();

    auto &container = layers[0];
    for (auto &cell: container.cells)
        cell.clear();

    container = layers[layers.size() - 1];
    for (auto &cell: container.cells)
        cell.clear();

}

void LinkedCell3D::preparePeriodic() {
    for (size_t i = 1; i < layers.size() - 1; ++i) {
        layers[i].updatePeriodic();
    }
    if (periodic.contains(Boundary::BACK)) {
        frontBackBoundary(-domain[2], layers.size() - 2, 0);
    }

    if (periodic.contains(Boundary::FRONT)) {
        frontBackBoundary(domain[2], 1, layers.size() - 1);
    }

}

void LinkedCell3D::frontBackBoundary(double to_add, size_t ind, size_t oth) {

    auto &lc = layers[ind];
    auto &counter = layers[oth];

    size_t i = LinkedCellContainer::mesh[0] + 1;
    for (; i < lc.cells.size() - mesh[0] - 1; ++i) {
        for (auto &p: lc[i]) {
            auto newP = p.getX();
            newP[2] += to_add;
            counter.simpleAdd(Particle(newP, p.getV(), p.getM(), p.getSigma(), p.getEpsilon(), p.getType()));

            mirrorHorizontal(newP, p, i, counter);
            mirrorVertical(newP, p, i, counter);
            mirrorDiagonal(newP, p, i, counter);
        }
    }


}

bool LinkedCell3D::side(size_t ind3D) {
    return (ind3D == 0 && periodic.contains(Boundary::FRONT)) ||
           (ind3D == layers.size() - 1 && periodic.contains(Boundary::BACK));

}

void LinkedCell3D::update(Particle &particle, size_t ind3D, size_t ind) {
    if (ind3D > 0 && ind3D < layers.size() - 1) {
        layers[ind3D].addParticle(particle);
        return;
    }

    if (side(ind3D) || layers[ind3D].side(ind)) {
        updatePeriodic(particle, ind3D, ind);
    }

}

void LinkedCell3D::updatePeriodic(Particle &p, size_t ind3D, size_t ind) {
    auto pos = p.getX();
    if (ind3D == 0) {
        pos[2] += domain[2];
        ind3D = layers.size() - 2;
    } else {
        pos[2] -= domain[2];
        ind3D = 1;
    }
    p.setX(pos);
    ind = layers[ind3D].mirror(p, ind);
    layers[ind3D][ind].addParticle(p);
}

void LinkedCell3D::setSize(double cutOff_arg, std::array<double, 3> &domain_arg) {
    domain = domain_arg;
    cutOff = cutOff_arg;
    for (size_t i = 0; i < 3; ++i) {
        mesh[i] = ceil(std::abs(domain_arg[i]) / cutOff_arg) + 2;
    }
    LinkedCellContainer::setDomain(domain_arg);
    LinkedCellContainer::setRCutOff(cutOff_arg);
    LinkedCellContainer::setMesh(mesh);

    layers.resize(mesh[2]);
    for (auto &layer: layers) {
        layer.setSize(cutOff_arg, domain_arg);
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

LinkedCellContainer &LinkedCell3D::operator[](size_t i) {
    return layers[i];
}

double LinkedCell3D::mirrorHorizontal(std::array<double, 3> &pos, Particle &p, size_t i, LinkedCellContainer &counter) {
    std::array<double, 3> newP = pos;
    if (LinkedCellContainer::leftBoundary(i)) {
        newP[0] += LinkedCellContainer::domain[0];
        counter.simpleAdd(Particle(newP, p.getV(), p.getM(), p.getSigma(), p.getEpsilon(), p.getType()));
        return domain[0];
    } else if (LinkedCellContainer::rightBoundary(i)) {
        newP[0] -= LinkedCellContainer::domain[0];
        counter.simpleAdd(Particle(newP, p.getV(), p.getM(), p.getSigma(), p.getEpsilon(), p.getType()));
        return -domain[0];
    }

    return 0.;
}

double LinkedCell3D::mirrorVertical(std::array<double, 3> &pos, Particle &p, size_t i, LinkedCellContainer &counter) {
    std::array<double, 3> newP = p.getX();
    if (LinkedCellContainer::bottomBoundary(i)) {
        newP[1] += LinkedCellContainer::domain[1];
        counter.simpleAdd(Particle(newP, p.getV(), p.getM(), p.getSigma(), p.getEpsilon(), p.getType()));
        return domain[1];
    } else if (LinkedCellContainer::topBoundary(i)) {
        newP[1] -= LinkedCellContainer::domain[1];
        counter.simpleAdd(Particle(newP, p.getV(), p.getM(), p.getSigma(), p.getEpsilon(), p.getType()));
        return -domain[1];
    }

    return 0;
}

void LinkedCell3D::mirrorDiagonal(std::array<double, 3> newP, Particle &p, size_t i, LinkedCellContainer &counter) {
    //mirror particles from corner cells to other corner cells
    if (LinkedCellContainer::topBoundary(i) && LinkedCellContainer::rightBoundary(i)) {
        //mirror boundary particle
        std::array<double, 3> to_add{-domain[0], -domain[1], 0};
        counter.simpleAdd(Particle(newP + to_add, p.getV(), p.getM(), p.getSigma(), p.getEpsilon(), p.getType()));

    } else if (LinkedCellContainer::topBoundary(i) && LinkedCellContainer::leftBoundary(i)) {
        //mirror boundary particle
        std::array<double, 3> to_add{domain[0], -domain[1], 0};
        counter.simpleAdd(Particle(newP + to_add, p.getV(), p.getM(), p.getSigma(), p.getEpsilon(), p.getType()));
    } else if (LinkedCellContainer::bottomBoundary(i) && LinkedCellContainer::leftBoundary(i)) {
        //mirror boundary particle
        std::array<double, 3> to_add{domain[0], domain[1], 0};
        counter.simpleAdd(Particle(p.getX() + to_add, p.getV(), p.getM(), p.getSigma(), p.getEpsilon(), p.getType()));

    } else if (LinkedCellContainer::bottomBoundary(i) && LinkedCellContainer::rightBoundary(i)) {
        //mirror boundary particle
        std::array<double, 3> to_add{-domain[0], domain[1], 0};
        counter.simpleAdd(Particle(newP + to_add, p.getV(), p.getM(), p.getSigma(), p.getEpsilon(), p.getType()));
    }
}


LinkedCell3D::LinkedCell3D() = default;
