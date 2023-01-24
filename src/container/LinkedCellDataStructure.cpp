//
// Created by lukas on 11.01.23.
//

#include "LinkedCellDataStructure.h"

LinkedCellDataStructure::~LinkedCellDataStructure() = default;
std::vector<Reflecting> LinkedCellDataStructure::conditions{};
std::set<Boundary> LinkedCellDataStructure::periodic{};
std::map<Boundary,Reflecting> LinkedCellDataStructure::frontBack{};

bool LinkedCellDataStructure::containsPeriodic(Boundary bound) {
    return periodic.find(bound) != periodic.end();
}

void LinkedCellDataStructure::addPeriodic(Boundary bound) {
    periodic.emplace(bound);
}

void LinkedCellDataStructure::addReflecting(Reflecting &&reflecting) {
    conditions.emplace_back(reflecting);
}



void LinkedCellDataStructure::clearBoundary() {
    periodic.clear();
    conditions.clear();
    frontBack.clear();
}

void LinkedCellDataStructure::addThird(Boundary bound, Reflecting &&ref) {
    frontBack.insert({bound, ref});
}
