//
// Created by lukas on 11.01.23.
//

#include "LinkedCellDataStructure.h"

LinkedCellDataStructure::~LinkedCellDataStructure() = default;

std::set<Boundary> LinkedCellDataStructure::periodic{};
std::map<Boundary,Reflecting> LinkedCellDataStructure::conditions{};

bool LinkedCellDataStructure::containsPeriodic(Boundary bound) {
    return periodic.find(bound) != periodic.end();
}

void LinkedCellDataStructure::addPeriodic(Boundary bound) {
    periodic.emplace(bound);
}

void LinkedCellDataStructure::addReflecting(Boundary bound,Reflecting &&reflecting) {
    conditions.insert({bound, reflecting});
}



void LinkedCellDataStructure::clearBoundary() {
    periodic.clear();
    conditions.clear();

}
