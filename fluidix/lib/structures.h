#ifndef STRUCTURES_H
#define STRUCTURES_H

#include "../lib/genome.h"
#include "../lib/nerveSystem.h"
#include <concurrent_queue.h>
#include <unordered_map>

struct Global {
    int nEggs;

    float dt; // integration time-step
    float pelletLifetime;
    int3 w;
    int nParticles;
    int nInitialOrganisms;
    int bufferSize;
    int nSteps;
    int3 initialOrganismDimensions;
    int nGenomeInputs; //Inputs x,y,z,d
    float interactionRange;
    int moveFactor;
    int repulsiveForce;
    float springForce; // spring constant
    float groundRepulsiveForce; // repulsive wall force
    float initialCellEnergy;
    float minCellEnergy;
    float minPelletEnergy;
    float energyParticleEnergy;
    float energyParticleRadius;
    float cellExistenceThreshold;
    float cellMetabolism;
    float cellDecayRate;
    float fluidDensity;
    float gravity;
} g;

struct Organism {
    Genome genome;
    NerveSystem nerveSystem;
    vector<int> cells;
    int parent;
};

enum CellType {
    Photo, Digest, Sting, Vascular, Fat, Sense, Egg,
    N_CELL_TYPES
};
enum ParticleType {
    Cell, Energy, Pellet, Buffer,
    N_PARTICLE_TYPES
};
enum Neigbour { Front, Right, Up, Back, Left, Down };

struct Particle {
    ParticleType particleType;
    xyz r, v, f;
    float color;
    float radius;
    float alpha;
    float density;
    float energy;
    float energyIn;
    float energyOut;
    float maxEnergy;
    float signal;
    int organism;
    bool toBuffer;
    int links[6];
    CellType type;
};

typedef std::unordered_map<int, Organism> OrganismMap;
typedef Concurrency::concurrent_queue<int> ParticleBuffer;

#endif