//
// Created by dominik on 29.10.22.
//
#pragma once

#include "../ParticleContainer.h"
#include "Particle.h"

namespace inputReader{
class InputTemplate {
public:

    virtual void readFile(ParticleContainer &particles, char *filename) = 0;
};
}
