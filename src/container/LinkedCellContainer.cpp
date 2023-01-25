//
// Created by lukas on 29.11.22.
//

#include <cmath>
#include <iostream>
#include "LinkedCellContainer.h"
#include "ParticleContainer.h"
#include "Reflecting.h"
#include "MolSimLogger.h"

std::array<double, 3> LinkedCellContainer::domain{};
double LinkedCellContainer::rcutoff{};
std::array<size_t, 3> LinkedCellContainer::mesh{};

void LinkedCellContainer::apply(std::function<void(Particle &)> fun) {
    for (size_t i = mesh[0] + 1; i < cells.size() - mesh[0] - 1; ++i) {
        cells[i].apply(fun);

        //skip halo
        if (i % mesh[0] == mesh[0] - 2)
            i += 2;
    }
}

void LinkedCellContainer::applyX(std::function<void(Particle &)> fun) {
    apply(fun);
    update();

}

void LinkedCellContainer::clearHalo() {
    for (auto &cell: halo)
        cell.get().clear();
}

void LinkedCellContainer::update() {
    size_t len = cells.size() - mesh[0] - 1;
    for (size_t i = mesh[0] + 1; i < len; ++i) {
        for (auto it = cells[i].begin(); it != cells[i].end();) {
            auto &p = *it;
            size_t ind = index(p);
            auto &pos = p.getX();

            //check that not completely outside
            if (pos[0] < -rcutoff || pos[0] >= domain[0] + rcutoff || pos[1] < -rcutoff ||
                pos[1] > domain[1] + rcutoff) {
                SPDLOG_LOGGER_INFO(MolSimLogger::logger(), "Particle at position ({}, {}, {}) removed", p.getX()[0],
                                   p.getX()[1], p.getX()[2]);
                it = cells[i].remove(it);
                continue;
            }
            //check if in the right cell
            if (ind == i) {
                ++it;
                continue;
            }
            //add to new cell
            update(p, ind);

            //remove from old cell
            it = cells[i].remove(it);
        }

        //skip halo
        if (i % mesh[0] == mesh[0] - 2)
            i += 2;
    }

    clearHalo();
}

size_t LinkedCellContainer::size() {

    size_t len = 0;
    for (size_t i = mesh[0] + 1; i < cells.size() - mesh[0] - 1; ++i) {
        len += cells[i].size();
        if (i % mesh[0] == mesh[0] - 2)
            i += 2;
    }
    return len;

}

void LinkedCellContainer::applyFBoundary(std::function<void(Particle &, Particle &)> &fun) {
    if(conditions.contains(Boundary::BOTTOM)) {
        auto & cond = conditions.at(Boundary::BOTTOM);
        for (auto &list: bottom_boundary) {
            for (auto &p: list.get()) {
                if (cond.check(p))
                    cond.apply(p, fun);
            }
        }
    }

    if(conditions.contains(Boundary::TOP)) {
        auto & cond = conditions.at(Boundary::TOP);
        for (auto &list: top_boundary) {
            for (auto &p: list.get()) {
                if (cond.check(p))
                    cond.apply(p, fun);
            }
        }
    }

    if(conditions.contains(Boundary::RIGHT)) {
        auto & cond = conditions.at(Boundary::RIGHT);
        for (auto &list: right_boundary) {
            for (auto &p: list.get()) {
                if (cond.check(p))
                    cond.apply(p, fun);
            }
        }
    }

    if(conditions.contains(Boundary::LEFT)) {
        auto & cond = conditions.at(Boundary::LEFT);
        for (auto &list: left_boundary) {
            for (auto &p: list.get()) {
                if (cond.check(p))
                    cond.apply(p, fun);
            }
        }
    }
}


void LinkedCellContainer::applyF(std::function<void(Particle &, Particle &)> fun) {
    updatePeriodic();
    size_t len = cells.size() - mesh[0] - 1;
#pragma omp parallel for schedule(dynamic, 1) num_threads(6) shared(cells)
    for (size_t i = 0; i < len; ++i) {
        auto &cell = cells[i];
        cell.applyF(fun);
        for (auto &p: cell) {
            auto partial = [&p, &fun](Particle &p2) { fun(p, p2); };

            //check if right neighbour exists
            rightNeighbour(i, partial);

            //check if upper neighbour exists
            upperNeighbour(i, partial);

            //check if upper right neighbour exists
            upperRightNeighbour(i, partial);

            //check is upper left neighbour exists
            upperLeftNeighbour(i, partial);
        }

    }

    applyFBoundary(fun);

}

size_t LinkedCellContainer::index(Particle &p) {
    auto &pos = p.getX();
    size_t x_ind, y_ind;

    if (pos[0] >= domain[0]) {
        x_ind = mesh[0] - 1;
    } else {
        x_ind = static_cast<int>(floor(pos[0] / rcutoff)) + 1;
    }

    if (pos[1] >= domain[1]) {
        y_ind = (mesh[1] - 1) * mesh[0];
    } else {
        y_ind = static_cast<int>(floor(pos[1] / rcutoff)) * mesh[0] + mesh[0];
    }
    return x_ind + y_ind;
}

void LinkedCellContainer::addParticle(Particle &&p) {
    Particle p1 = p;
    size_t ind = index(p1);
    auto &pos = p.getX();
    if (0 <= pos[0] && pos[0] < domain[0] && 0 <= pos[1] && pos[1] < domain[1]) {
        cells[ind].addParticle(p);

    }
}

std::vector<ParticleContainer> &LinkedCellContainer::getCells() {
    return cells;
}

void LinkedCellContainer::setUpRef() {
    size_t i = 0;
    size_t len = cells.size();

    //add first row of halo cells
    for (; i < mesh[0]; ++i)
        halo.emplace_back(std::ref(cells[i]));

    //additional halo cell on the left sie
    halo.emplace_back((std::ref(cells[i])));
    ++i;

    left_boundary.emplace_back(std::ref(cells[i]));
    //add all boundary cell in the bottom boundary
    size_t j = i;
    for (; j < i + mesh[0] - 2; ++j)
        bottom_boundary.emplace_back(std::ref(cells[j]));
    i = j;

    right_boundary.emplace_back(std::ref(cells[i - 1]));


    //add halo cell at the right border
    halo.emplace_back(std::ref(cells[i]));
    ++i;

    //add right and left halo cell as well as right an left boundary cell
    //or nothing if mesh[1] = 1
    for (; i < len - 2 * mesh[0]; ++i) {
        //right or left halo cell
        if (i % mesh[0] == 0 || i % mesh[0] == mesh[0] - 1) {
            halo.emplace_back(std::ref((cells[i])));
            continue;
        }

        //right or left boundary cell
        //if mesh[0] only middle cell
        if (i % mesh[0] == 1) {
            left_boundary.emplace_back(std::ref(cells[i]));
        }else if (i % mesh[0] == mesh[0] - 2){
            right_boundary.emplace_back(std::ref(cells[i]));
        }

    }

    //add halo cell at left boundary from last boundary row
    //or left halo cell of uppest row if mesh[1] = 1
    halo.emplace_back(std::ref(cells[i]));
    ++i;

    left_boundary.emplace_back(std::ref(cells[i]));
    //add upper boundary cells or nothing if mesh[1] = 1
    j = i;
    for (; j < len - mesh[0] - 1; ++j)
        top_boundary.emplace_back(std::ref(cells[j]));
    i = j;

    right_boundary.emplace_back(std::ref(cells[i - 1]));
    //add halo cell at right border or halo cell in the middle if mesh[1] = 1
    halo.emplace_back(std::ref(cells[i]));
    ++i;

    //add halo cells at upper border
    for (; i < len; ++i)
        halo.emplace_back(std::ref(cells[i]));
}

void LinkedCellContainer::setRCutOff(double rcutoff_arg) {
    rcutoff = rcutoff_arg;
}

void LinkedCellContainer::setDomain(std::array<double, 3> &domain_arg) {
    domain = domain_arg;
}

void LinkedCellContainer::setSize(double rcutoff_arg, std::array<double, 3> &domain_arg) {
    setRCutOff(rcutoff_arg);
    setDomain(domain_arg);
    for (size_t i = 0; i < 2; ++i) {
        mesh[i] = static_cast<size_t>(ceil(std::abs(domain_arg[i]) / rcutoff_arg)) + 2;
    }
    setUp();
}


const std::vector<std::reference_wrapper<ParticleContainer>> &LinkedCellContainer::getHalo() const {
    return halo;
}

std::vector<std::reference_wrapper<ParticleContainer>> LinkedCellContainer::getBoundary() const {
    std::vector<std::reference_wrapper<ParticleContainer>> bound;
    bound.insert(bound.cend(),bottom_boundary.begin(), bottom_boundary.end());
    bound.insert( bound.cend(), top_boundary.begin(), top_boundary.end());
    bound.insert(bound.cend(), right_boundary.begin(), right_boundary.end());
    bound.insert(bound.cend(),left_boundary.begin(), left_boundary.end());
    return bound;
}

std::array<double, 3> &LinkedCellContainer::getDomain() {
    return domain;
}

LinkedCellContainer::LinkedCellContainer() = default;

LinkedCellContainer::~LinkedCellContainer() = default;

void LinkedCellContainer::rightNeighbour(size_t i, const std::function<void(Particle &)> &partial) {
    if (mesh[0] <= 1)
        return;

    if ((i + 1) % mesh[0] > 0) {
        auto &neighbour = cells[i + 1];
        neighbour.apply(partial);
    }
}

void LinkedCellContainer::upperNeighbour(size_t i, const std::function<void(Particle &)> &partial) {
    if (mesh[1] <= 1)
        return;

    if (i < cells.size() - mesh[0] - 1) {

        auto &neighbour = cells[i + mesh[0]];
        neighbour.apply(partial);
    }
}

void LinkedCellContainer::upperLeftNeighbour(size_t i, const std::function<void(Particle &)> &partial) {
    if (mesh[0] <= 1 || mesh[1] <= 1)
        return;

    if (i < cells.size() - mesh[0] - 1 && i % mesh[0] > 0) {

        auto &neighbour = cells[i + mesh[0] - 1];
        neighbour.apply(partial);
    }
}

void LinkedCellContainer::upperRightNeighbour(size_t i, const std::function<void(Particle &)> &partial) {
    if (mesh[0] <= 1 || mesh[1] <= 1)
        return;

    if (i < cells.size() - mesh[0] - 1 && (i + 1) % mesh[0] > 0) {
        auto &neighbour = cells[i + mesh[0] + 1];
        neighbour.apply(partial);
    }
}


bool LinkedCellContainer::side(size_t ind) {
    return (ind % mesh[0] == 0 && containsPeriodic(Boundary::LEFT)) ||
           (ind % mesh[0] == mesh[0] - 1 && containsPeriodic(Boundary::RIGHT)) ||
           (ind < mesh[0] && containsPeriodic((Boundary::BOTTOM))) ||
           (ind >= cells.size() - mesh[0] && containsPeriodic(Boundary::TOP));
}

size_t LinkedCellContainer::mirror(Particle &p, size_t ind) {

    //move to right boundary
    if (ind % mesh[0] == 0) {
        std::array<double, 3> to_add{};
        to_add[0] = domain[0];
        p.setX(p.getX() + to_add);
        ind = index(p);
    }

    //move to right boundary
    if (ind % mesh[0] == mesh[0] - 1) {
        std::array<double, 3> to_add{};
        to_add[0] = -domain[0];
        p.setX(p.getX() + to_add);
        ind = index(p);
    }

    //move to upper boundary
    if (ind < mesh[0]) {
        std::array<double, 3> to_add{};
        to_add[1] = domain[1];
        p.setX(p.getX() + to_add);
        ind = index(p);
    }

    //move to lower boundary
    if (ind >= cells.size() - mesh[0]) {
        std::array<double, 3> to_add{};
        to_add[1] = -domain[1];
        p.setX(p.getX() + to_add);
        ind = index(p);
    }

    return ind;
}

void LinkedCellContainer::update(Particle &p, size_t ind) {
    if (side(ind)) {
        ind = mirror(p, ind);
    }
    cells[ind].addParticle(p);

}

void LinkedCellContainer::mirrorPeriodic() {

    //mirror particles from corner cells to other corner cells
    if (containsPeriodic(Boundary::BOTTOM) || containsPeriodic(Boundary::LEFT)) {
        //mirror boundary particle
        std::array<double, 3> to_add{-domain[0], -domain[1], 0};
        auto & cell = top_boundary[top_boundary.size() - 1];
        for(auto & p: cell.get()) {
            simpleAdd(Particle(p.getX() + to_add, p.getV(), p.getM(), p.getSigma(), p.getEpsilon(), p.getType()));
        }

    }
    if (containsPeriodic(Boundary::BOTTOM) || containsPeriodic(Boundary::RIGHT)) {
        //mirror boundary particle
        std::array<double, 3> to_add{domain[0], -domain[1], 0};
        auto & cell = top_boundary[0];
        for(auto & p: cell.get()) {
            simpleAdd(Particle(p.getX() + to_add, p.getV(), p.getM(), p.getSigma(), p.getEpsilon(), p.getType()));
        }
    }
    if (containsPeriodic(Boundary::TOP) || containsPeriodic(Boundary::RIGHT)) {
        //mirror boundary particle
        std::array<double, 3> to_add{domain[0], domain[1], 0};
        auto & cell = bottom_boundary[0];
        for(auto & p: cell.get()) {
            simpleAdd(Particle(p.getX() + to_add, p.getV(), p.getM(), p.getSigma(), p.getEpsilon(), p.getType()));
        }
    }
    if (containsPeriodic(Boundary::TOP) || containsPeriodic(Boundary::LEFT)) {
        //mirror boundary particle
        std::array<double, 3> to_add{-domain[0], domain[1], 0};
        auto & cell = right_boundary[0];
        for(auto & p: cell.get()) {
            simpleAdd(Particle(p.getX() + to_add, p.getV(), p.getM(), p.getSigma(), p.getEpsilon(), p.getType()));
        }
    }
}

void LinkedCellContainer::simpleAdd(Particle &&p) {
    size_t ind = index(p);
    cells[ind].addParticle(p);

}

bool LinkedCellContainer::bottomBoundary(size_t ind) {
    return ind > mesh[0] && ind < 2 * mesh[0] - 1;
}

bool LinkedCellContainer::leftBoundary(size_t ind) {
    return ind % mesh[0] == 1;
}

bool LinkedCellContainer::rightBoundary(size_t ind) {
    return ind % mesh[0] == mesh[0] - 2;
}

bool LinkedCellContainer::topBoundary(size_t ind) {
    size_t len = mesh[0] * mesh[1] - mesh[0];
    return ind < len - 1 && ind > len - mesh[0];
}

void LinkedCellContainer::updatePeriodic() {
    //mirror particles for periodic boundary
    if(periodic.contains(Boundary::RIGHT)) {
        for (auto &cell: left_boundary) {
            mirrorBoundary(cell, {domain[0], 0, 0});
        }
    }

    if(periodic.contains(Boundary::LEFT)) {
        for (auto &cell: right_boundary) {
            mirrorBoundary(cell, {-domain[0], 0, 0});
        }
    }

    if(periodic.contains(Boundary::TOP)) {
        for (auto &cell: bottom_boundary) {
            mirrorBoundary(cell, {0, domain[1], 0});
        }
    }

    if(periodic.contains(Boundary::BOTTOM)) {
        for (auto &cell: top_boundary) {
            mirrorBoundary(cell, {0, -domain[1], 0});
        }
    }

    mirrorPeriodic();

}

void LinkedCellContainer::addParticle(Particle &p) {
    size_t ind = index(p);
    auto &pos = p.getX();
    if (0 <= pos[0] && pos[0] < domain[0] && 0 <= pos[1] && pos[1] < domain[1]) {
        cells[ind].addParticle(p);

    }
}

ParticleContainer &LinkedCellContainer::operator[](size_t i) {
    return cells[i];
}

void LinkedCellContainer::forceTwoD(ParticleContainer &particles, size_t ind,
                                    std::function<void(Particle &, Particle &)> fun) {
    particles.applyF(fun);
    for (auto &p: particles) {
        auto partial = [&p, &fun](Particle &p2) { fun(p, p2); };
        rightNeighbour(ind, partial);
        upperLeftNeighbour(ind, partial);
        upperRightNeighbour(ind, partial);
        upperNeighbour(ind, partial);
    }

}

void LinkedCellContainer::setUp() {
    size_t len = 1;
    for (size_t i = 0; i < 2; ++i) {
        len *= mesh[i] == 0 ? 1 : mesh[i];
    }
    cells.resize(len);
    setUpRef();
}

void LinkedCellContainer::setMesh(std::array<size_t, 3> &mesh_arg) {
    mesh = mesh_arg;
}

void LinkedCellContainer::clearBoundary() {
    periodic.clear();
    conditions.clear();
}

void LinkedCellContainer::applyPar(std::function<void(Particle &)> fun) {
    size_t len = mesh[0] * mesh[1];

#pragma omp parallel for schedule(dynamic , 1) num_threads(6)
    for (size_t i = mesh[0] + 1; i < len - mesh[0] - 1; ++i) {
        if (i % mesh[0] == 0 || i % mesh[0] == mesh[0] - 1)
            continue;
        cells[i].apply(fun);
    }
}

void LinkedCellContainer::mirrorBoundary(ParticleContainer &par, std::array<double, 3> &&to_add) {
    for(auto& p : par){
        simpleAdd(Particle(p.getX() + to_add, p.getV(), p.getM(), p.getSigma(), p.getEpsilon(), p.getType()));
    }
}

