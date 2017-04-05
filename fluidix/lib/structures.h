#ifndef STRUCTURES_H
#define STRUCTURES_H

#include "../lib/genome.h"
#include "../lib/nerveSystem.h"
//#include <concurrent_queue.h>
#include <queue>
#include <unordered_map>

struct Global {
    float dt; // integration time-step
    int3 w;
    int nParticles;
    int nInitialOrganisms;
    int bufferSize;
    int energyParticleCount;
    int nSteps;
    int nGenomeInputs; //Inputs x,y,z,d
    float interactionRange;
    int moveFactor;
    float moveCost;
    float genomeCost;
    float nerveCost;
    int repulsiveForce;
    float springForce; // spring constant
    float groundRepulsiveForce; // repulsive wall force
    float initialCellEnergy;
    float minCellEnergy;
    float minDetritusEnergy;
    float energyParticleEnergy;
    float energyParticleRadius;
    float cellExistenceThreshold;
    float cellMetabolism;
    float cellDecayRate;
    float fluidDensity;
    float gravity;

   int nDetritus;
   int nBuffer;
   int nEnergy;
   int nCells;
} g;

struct Organism {
    Genome genome;
    NerveSystem nerveSystem;
    vector<int> cells;
    int parent;
    float health;
};

enum CellType {
    Photo, Digest, Sting, Vascular, Fat, Sense, Egg,
    N_CELL_TYPES
};
enum ParticleType {
    Cell, Energy, Detritus, Buffer,
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
    float metabolism;
    int organism;
    bool toBuffer;
    int links[6];
    CellType type;
};

typedef std::unordered_map<int, Organism> OrganismMap;
//typedef Concurrency::concurrent_queue<int> ParticleBuffer;
typedef std::queue<int> ParticleBuffer;

#endif