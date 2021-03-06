//#include "fluidix.h"
#define PARTICLE_BITS 28
#include "C:\Program Files (x86)\Fluidix\include\fluidix.h"
#include "../lib/genome.h"
#include <queue>
#include <ppl.h>
#include <concurrent_queue.h>
#include <unordered_map>

#define DT 0.01f // integration time-step
#define PELLET_LIFETIME 5.0f

#define W make_int3(200, 100, 200)
#define N 100000
#define N_ORIGIN_CELLS 1000
#define N_INITIAL_BUFFER 50000
#define MAX_NERVE_CONNECTIONS 100000
#define N_STEPS 10000000
#define INITIAL_ORGANISM_DIMENSIONS make_int3(1, 1, 1)

//Inputs x,y,z,d:
#define N_INPUTS 4

#define RANGE 10.0f
#define MOVE_FACTOR 10

#define REPULSION_FORCE 300
#define SPRING_K 150.0f // spring constant
#define WALL 100.0f // repulsive wall force

#define CELL_INITIAL_ENERGY 2.0f
#define CELL_MIN_ENERGY 1.0f
#define PELLET_MIN_ENERGY 0.01f
#define DIVISION_ENERGY (CELL_MIN_ENERGY * 2.5f)
#define ENERGY_PARTICLE_ENERGY 0.5f
#define MAX_CELL_DIVISIONS 100
#define CELL_EXISTENCE_THRESHOLD 0.0f
#define NERVE_LEARNING_RATE 0.1f

#define CELL_METABOLISM 0.02f
#define CELL_DECAY_RATE 0.01f

#define FLUID_DENSITY 1.0f
#define G -9.81f

#define GREEN   0.5f
#define RED     1.0f
#define YELLOW  0.7f
#define BLUE    0.0f
#define CYAN    0.3f
#define ORANGE  0.8f

#define transmitFloat(a, b, f) {addFloat(a, -f); addFloat(b, f);}
#define isWeirdParticle(p) (p.r.x != p.r.x || p.r.y != p.r.y || p.r.z != p.r.z)
#define isWeirdParticlePointer(p) (p->r.x != p->r.x || p->r.y != p->r.y || p->r.z != p->r.z)

using namespace std;
using namespace concurrency;

typedef std::unordered_multimap<int, int> intMap;

enum CellType {
    Photo, Digest, Sting, Vascular, Fat, Sense, Motor, Ballast, Egg, Neuron, 
    N_CELL_TYPES
};
enum ParticleType {
    Cell, Energy, Pellet, Buffer,
    N_PARTICLE_TYPES
};
enum Neigbour { Front, Right, Up, Back, Left, Down };

int currGenomeIndex;

struct Global {
    int nCells;
    float *nerveWeights;
} g;

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
    int nDivisions;
    float preSignal;
    float postSignal;
    int organism;
    bool toBuffer;
    bool toReproduce;
    Particle *origin;
    int links[6];
    bool missingNeighbour[6];
    CellType type;
};

#define turnIntoPellet(p) {         \
    p.particleType = Pellet;        \
    p.energy = PELLET_LIFETIME;     \
    p.density = FLUID_DENSITY * 2;  \
    p.alpha = 0.5f;                 \
    p.organism = -1;                \
}                                   \

#define turnIntoEnergy(p) {         \
    p.particleType = Energy;        \
    p.r = make_xyz(                 \
        rnd_uniform()*W.x,          \
        rnd_uniform()*W.y + W.y,    \
        rnd_uniform()*W.z           \
    );                              \
    p.color = 0.7f;                 \
    p.energy = ENERGY_PARTICLE_ENERGY; \
    p.preSignal = 0.0f;             \
    p.postSignal = 0.0f;            \
    p.alpha = 0.3f;                 \
    p.radius = 0.5f;                \
    p.density = 10.0f; 	  				\
}

#define turnIntoBuffer(p) {         \
    p.particleType = Buffer;        \
    p.r = make_xyz(                 \
        rnd_uniform()*W.x,          \
        rnd_uniform()*W.y - W.y,    \
        rnd_uniform()*W.z           \
    );                              \
    p.density = FLUID_DENSITY * 2;  \
    p.alpha = 0.5f;                 \
    p.color = 0.5f;                 \
    p.radius = 1.0f;                \
    p.organism = -1;                \
}

FUNC_EACH(init,
    turnIntoEnergy(p);
)

FUNC_GENERIC(initNerveWeights,
    g.nerveWeights[i] = 1.0f;
)

FUNC_EACH(integrate,
    p.v += p.f * DT; //Mass?
    p.r += p.v * DT;
    p.f = make_xyz(0, 0, 0);
    p.v *= 0.97f;
)

FUNC_EACH(handleEnergy,
    switch (p.particleType) {
    case Cell:
        p.energy -= CELL_METABOLISM * DT;
        if (p.energy < CELL_MIN_ENERGY )
            turnIntoPellet(p)
        else if (p.energy > p.maxEnergy)
            p.energy -= (p.energy - p.maxEnergy) * 0.1f;
        break;
    case Pellet:
        p.energy -= CELL_DECAY_RATE * DT;
        if (p.energy <= PELLET_MIN_ENERGY)
            p.toBuffer = true;
        break;
    }
    //p.color = mapf(p.energy, 0.0f, 3.0f, 0.0f, 1.0f);
)

FUNC_EACH(buoyancy,
    //if (p.particleType == Cell && p.type == Ballast)
    //    p.density = clamp(p.density + p.signal, 0.5f, 2.0f);
    float volume = p.radius * p.radius * PI;
    p.f.y += (p.density - FLUID_DENSITY) * G * volume;
)

FUNC_EACH(moveParticle,
    if (p.particleType == Cell && p.type == Motor) {
        xyz f = make_xyz(0, 0, 0);
        if (p.links[Left]  < 0)  f.x += 1;
        if (p.links[Up]    < 0)  f.y += 1;
        if (p.links[Back]  < 0)  f.z += 1;
        if (p.links[Right] < 0)  f.x -= 1;
        if (p.links[Down]  < 0)  f.y -= 1;
        if (p.links[Front] < 0)  f.z -= 1;
        f = xyz_norm(f) * MOVE_FACTOR; //* rnd_normal(); //clamp(p.preSignal, -1.0f, 1.0f);

        addVector(p.f, f);
    }
)

FUNC_EACH(reproduction,
    if (p.particleType == Cell && p.type == Egg)
        if (p.energy >= p.maxEnergy){
            p.toReproduce = true;
            printf("toReproduce!!!\n");
        }
)

FUNC_EACH(resetMissingNeighbours,
    for (int i = 0; i < 6; i++)
        p.missingNeighbour[i] = true;
)

FUNC_EACH(handleMissingNeighbours,
    for (int i = 0; i < 6; i++)
    if (p.missingNeighbour[i] && p.links[i] != -1){
            p.toBuffer = true;
            break;
        }
)

#define sigm(x) x / (1.0f + abs(x))
FUNC_LINK(transmitNerveSignals,
    float w = g.nerveWeights[i];
    if ((p1.type == Neuron && p2.type == Sense) ||
        (p2.type == Neuron && p1.type == Motor) ||
        (p1.type == Neuron && p1.type == Neuron)
        )
    {
        float signal = w * p1.preSignal;
        addFloat(p2.postSignal, signal);
        //addFloat(g.nerveWeights[i], signal);
    }
    if ((p2.type == Neuron && p1.type == Sense) ||
        (p1.type == Neuron && p2.type == Motor) ||
        (p2.type == Neuron && p2.type == Neuron)
        )
    {
        float signal = w * p1.preSignal;
        addFloat(p1.postSignal, signal);
        //addFloat(g.nerveWeights[i], signal);
    }
)

FUNC_EACH(nerveActivationFunction,
    p.preSignal = sigm(p.postSignal);
    //p.color = p.preSignal;
    //if ((p.type == Neuron || p.type == Sense || p.type == Motor))
    //    printf("Signal of type %i is equal to: %.2f\n", p.type, p.preSignal);
    p.postSignal = 0.0f;
)

FUNC_LINK(updateNerveWeights,
    float w = g.nerveWeights[i];
    //if (w !=1.0) printf("Weight %i is %.2f\n", i, w);

    if ((p1.type == Neuron && p2.type == Sense) ||
        (p2.type == Neuron && p1.type == Motor) ||
        (p1.type == Neuron && p1.type == Neuron)
        )
    {
        float x = p1.preSignal;
        float y = p2.preSignal;

        float dw = NERVE_LEARNING_RATE * (x*y - y*y*w);

        //Oja learning rule:
        addFloat(g.nerveWeights[i], dw);
    }
    if ((p2.type == Neuron && p1.type == Sense) ||
        (p1.type == Neuron && p2.type == Motor) ||
        (p2.type == Neuron && p2.type == Neuron)
        )
    {
        float x = p2.preSignal;
        float y = p1.preSignal;

        float dw = NERVE_LEARNING_RATE * (x*y - y*y*w);

        //Oja learning rule:
        addFloat(g.nerveWeights[i], dw);
    }
)

// bouncing hard wall boundary condition
FUNC_EACH(boundary,
    if (p.particleType == Cell)
    addInteger(g.nCells, 1);

    // Check for wierd 1.#R values... NaN?
    if (isWeirdParticle(p)) {
        //printf("Weird particle!!! type (%i:%i)\n", p.particleType, p.type);
        turnIntoBuffer(p);
        p.toBuffer = true;
    } if (p.particleType != Buffer) {
        if (p.r.x < 0)   p.f.x += WALL * (0 - p.r.x);
        if (p.r.x > W.x) p.f.x += WALL * (W.x - p.r.x);
        if (p.r.z < 0)   p.f.z += WALL * (0 - p.r.z);
        if (p.r.z > W.x) p.f.z += WALL * (W.x - p.r.z);

        if (p.particleType == Energy){
            if (p.r.y < 0) {
                p.toBuffer = true;
            }
        }
        else {
            if (p.r.y < 0)   { p.v.y = 0.9f * (0 - p.r.y)   / DT; p.r.y = 0; }
            if (p.r.y > W.y) { p.v.y = 0.9f * (W.y - p.r.y) / DT; p.r.y = W.y; }
        }
    }
    else if (p.r.y < -2 * W.y) p.r.y += W.y;
)

#define consumeParticle(a, b) {             \
    addFloat(a.energy, b.energy);    \
    b.energy = 0;                           \
    b.toBuffer = true;                      \
}

FUNC_PAIR(particlePair,
    if (p1.particleType != Buffer && p2.particleType != Buffer) {
        xyz f = u * maxf(
            (REPULSION_FORCE * (1 - dr / (p1.radius + p2.radius))),
            0
        );
        if (p1.particleType == Cell &&
            p2.particleType == Cell)
        {
            //Cells from the same organism
            if (p1.organism == p2.organism) {
                bool neighbours = false;
                for (int n = 0; n < 6; n++) {
                    if (p1.links[n] == p2_index) {
                        p1.missingNeighbour[n] = false;
                        neighbours = true;
                    }
                    if (p2.links[n] == p1_index) {
                        p1.missingNeighbour[n] = false;
                        neighbours = true;
                    }
                }
                if (neighbours) {
                    //Spring force between neighbours
                    f = -u * ((dr - (p1.radius + p2.radius)/2) * SPRING_K);

                    //Signalling between cells of same organism
                    //float meanSignal = (p1.signal + p2.signal) / 2;
                    //p1.signal = p2.signal = meanSignal;

                    //Energy transmission
                    float p1Surplus = maxf(p1.energy - CELL_MIN_ENERGY, 0);
                    float p2Surplus = maxf(p2.energy - CELL_MIN_ENERGY, 0);
                    transmitFloat(p1.energy, p2.energy, p1Surplus * p1.energyOut * p2.energyIn);
                    transmitFloat(p2.energy, p1.energy, p2Surplus * p2.energyOut * p2.energyIn);
                }
            }
            //Cells from different organisms
            else {
                //Kill the other cell if you are sting
                if (p1.type == Sting && dr <= (p1.radius + p2.radius)) {
                    //turnIntoPellet(p2);
                    transmitFloat(p2.energy, p1.energy, 0.01f);
                }
                if (p2.type == Sting && dr <= (p1.radius + p2.radius)) {
                    //turnIntoPellet(p2);
                    transmitFloat(p1.energy, p2.energy, 0.01f);
                }
            }
        }
        //If p1 is a cell
        else if (p1.particleType == Cell) {
            if ((p1.type == Photo && p2.particleType == Energy) ||
                (p1.type == Digest  && p2.particleType == Pellet)
                ) consumeParticle(p1, p2)
        }
        //If p2 is a cell
        else if (p2.particleType == Cell) {
            if ((p2.type == Photo && p1.particleType == Energy) ||
                (p2.type == Digest  && p1.particleType == Pellet)
                ) consumeParticle(p2, p1)
        }

        if (p1.particleType == Cell && p1.type == Sense) addFloat(p1.preSignal, 1.0f);
        if (p2.particleType == Cell && p2.type == Sense) addFloat(p1.preSignal, 1.0f);

        addVector(p1.f, f);
        addVector(p2.f, -f);

        //p1.color = p1.signal;
        //p2.color = p2.signal;
    }
)

void setDefaultCellValues(Particle *cell) {
    cell->alpha = 1.0f;
    cell->radius = 1.0f;
    cell->energy = CELL_INITIAL_ENERGY;
    cell->density = FLUID_DENSITY * 1.10f;
    cell->particleType = Cell;
}

bool applyPhenotype(vector<float> output, Particle *cell) {
    // If cell should not exist, return
    if (output[N_CELL_TYPES + 1] < CELL_EXISTENCE_THRESHOLD)
        return false;

    float max = output[0]; cell->type = (CellType)0;
    for (int j = 1; j<N_CELL_TYPES; j++) {
        if (output[j] > max) {
            max = output[j];
            cell->type = (CellType)j;
        }
    }
    switch (cell->type) {
    case Photo:
        cell->color = GREEN;
        cell->energyIn = 0.01f;
        cell->energyOut = 0.6f;
        cell->maxEnergy = 3.0f;
        break;
    case Digest:
        cell->color = RED;
        cell->energyIn = 0.01f;
        cell->energyOut = 0.6f;
        cell->maxEnergy = 3.0f;
        break;
    case Fat:
        cell->color = YELLOW;
        cell->energyIn = 1.0f;
        cell->energyOut = 0.01f;
        cell->maxEnergy = 10.0f;
        break;
    case Motor:
        cell->color = 0.4f;
        cell->energyIn = 1.0f;
        cell->energyOut = 0.0f;
        cell->maxEnergy = 3.0f;
        break;
    case Sense:
        cell->color = BLUE;
        cell->energyIn = 1.0f;
        cell->energyOut = 0.0f;
        cell->maxEnergy = 3.0f;
        break;
    case Neuron:
        cell->color = nanf("gray");
        cell->energyIn = 1.0f;
        cell->energyOut = 0.0f;
        cell->maxEnergy = 3.0f;
        break;
    case Ballast:
        cell->color = CYAN;
        cell->energyIn = 1.0f;
        cell->energyOut = 0.0f;
        cell->maxEnergy = 3.0f;
        break;
    case Egg:
        cell->color = ORANGE;
        cell->energyIn = 1.0f;
        cell->energyOut = 0.0f;
        cell->maxEnergy = 1000.0f;
        break;
    case Vascular:
        cell->color = 0.2;
        cell->energyIn = 1.0f;
        cell->energyOut = 0.2f;
        cell->maxEnergy = 1.0f;
        break;
    case Sting:
        cell->color = nanf("");
        cell->energyIn = 1.0f;
        cell->energyOut = 0.0f;
        cell->maxEnergy = 3.0f;
        break;
    default:
        cell->energyIn = 1.0f;
        cell->energyOut = 0.0f;
        cell->maxEnergy = 3.0f;
    }
    if (cell->origin == cell) {
        cell->energyIn = 1.0f;
        //cell->energyOut = 0.0f;
        //cell->maxEnergy = 3.0f;
    }
    cell->nDivisions = output[N_CELL_TYPES] * MAX_CELL_DIVISIONS;
    return true;
}

void disconnectCell(Particle *p, int cell, int code) {
    for (int i = 0; i < 6; i++) {
        if (p[cell].links[i] >= 0)
            p[p[cell].links[i]].links[(i + 3) % 6] = code;
    }
    turnIntoBuffer(p[cell]);
}
void emptyCellPos(Particle *p, int cell) {
    disconnectCell(p, cell, -1);
}
void deadCellPos(Particle *p, int cell) {
    disconnectCell(p, cell, -2);
}

void addNerve(intMap *nerveMap, Fluidix<> *fx, int linkSet, int particleSet, int iFrom, int iTo)
{
    fx->addLink(linkSet, particleSet, iFrom, particleSet, iTo);
    nerveMap->insert(make_pair(iFrom, iTo));
    nerveMap->insert(make_pair(iTo, iFrom));
}

int getIdxFromCoord(int x, int y, int z, int3 br)
{
    x += br.x; y += br.y; z += br.z;
    //int lX = 2 * br.z + 1;
    int lY = 2 * br.y + 1;
    int lZ = 2 * br.z + 1;
    return x*lY*lZ + y*lZ + z;
}
#define iFromCoord(x,y,z) cellBuff.at(getIdxFromCoord(x,y,z,br))

// Initialize new organism
bool spawnOrganism(
    Fluidix<> *fx, int particleSet, int linkSet,
    intMap *nerveMap,
    xyz origin, concurrent_queue<int> *particleBuffer, 
    Particle *p, Genome *genome, Genome *allGenomes)
{
    int nParticlesNeeded = genome->getMaxCellsReq();
    if (nParticlesNeeded > particleBuffer->unsafe_size()) {
        printf("Not enought particles in buffer\n");
        return false;
    }
    vector<int> cellBuff;
    while (nParticlesNeeded) {
        int particle;
        if (particleBuffer->try_pop(particle)) {
            cellBuff.push_back(particle);
            nParticlesNeeded--;
        }
    }

    int3 br = genome->getBoundingRadius();

    int organism = currGenomeIndex++;
    vector<int> removedCells;
    vector<int> neurons;
    vector<int> allNerveCells;

    for (int x = -br.x; x <= br.x; x++)
    for (int y = -br.y; y <= br.y; y++)
    for (int z = -br.z; z <= br.z; z++) {
        Particle *cell = &p[iFromCoord(x, y, z)];
        cell->organism = organism;
        cell->r = origin + make_xyz(x, y, z);
        cell->toReproduce = false;
        cell->energy = CELL_INITIAL_ENERGY;
        setDefaultCellValues(cell);

        vector<float> input;
        input.push_back(x);
        input.push_back(y);
        input.push_back(z);
        input.push_back(xyz_len(make_xyz(x,y,z)));

        vector<float> output = genome->getOutput(input);

        if (applyPhenotype(output, cell)) {
            Genome gNew = Genome(*genome);
            gNew.mutate();
            allGenomes[iFromCoord(x, y, z)] = gNew;
            cell->links[Left] = x + 1 < br.x ? iFromCoord(x + 1, y, z) : -1;
            cell->links[Up] = y + 1 < br.y ? iFromCoord(x, y + 1, z) : -1;
            cell->links[Back] = z + 1 < br.z ? iFromCoord(x, y, z + 1) : -1;
            cell->links[Right] = x - 1 >= 0 ? iFromCoord(x - 1, y, z) : -1;
            cell->links[Down] = y - 1 >= 0 ? iFromCoord(x, y - 1, z) : -1;
            cell->links[Front] = z - 1 >= 0 ? iFromCoord(x, y, z - 1) : -1;

            if (cell->type == Neuron)
                neurons.push_back(iFromCoord(x, y, z));
            if (cell->type == Neuron ||
                cell->type == Sense  || 
                cell->type == Motor)
               allNerveCells.push_back(iFromCoord(x, y, z));
        }
        else
            removedCells.push_back(iFromCoord(x, y, z));
    }
    for (int i : removedCells)
        emptyCellPos(p, i);

    for (int iFrom : neurons)
        for (int iTo : allNerveCells)
            addNerve(nerveMap, fx, linkSet, particleSet, iFrom, iTo);

    return true;
}

#define printP(chr, p, i) printf("%c\tp[%i].r=(%.2f, %.2f, %.2f)\n", chr, i, p.r.x, p.r.y, p.r.z)

void generateTerrain(Fluidix<> *fx){
    normal_distribution<float> rndNormal(0.0f, W.y * 0.1);

    int meshParticles   = fx->createParticleSet(W.x * W.z);
    int meshLinks       = fx->createLinkSet();
    Particle *mesh      = fx->getParticleArray(meshParticles);

    for (int x = 0; x < W.x; x++)
    for (int z = 0; z < W.z; z++){
        int i = x*W.z + z;
        float y = rndNormal(rndGen);
        mesh[i].r = make_xyz(x, y, z);

        int n = (x+1)*W.z + z;
        int s = (x-1)*W.z + z;
        int e = x*W.z + (z+1);
        int w = x*W.z + (z-1);
        fx->addLink(meshLinks, meshParticles, i, meshParticles, n);
        fx->addLink(meshLinks, meshParticles, i, meshParticles, s);
        fx->addLink(meshLinks, meshParticles, i, meshParticles, e);
        fx->addLink(meshLinks, meshParticles, i, meshParticles, w);
    }
}

int main() {
    Fluidix<> *fx = new Fluidix<>(&g);
    int pSet = fx->createParticleSet(N);

    //generateTerrain(fx);

    fx->createGlobalArray(&g.nerveWeights, MAX_NERVE_CONNECTIONS);
    fx->runGeneric(initNerveWeights(), MAX_NERVE_CONNECTIONS);

    currGenomeIndex = 0;
    g.nCells = 0;
    Genome *genomes = new Genome[N];
    
    int lSet = fx->createLinkSet();
    intMap nerveMap;

    fx->runEach(init(), pSet);
    concurrent_queue<int> particleBuffer;
    Particle *p = fx->getParticleArray(pSet);

    for (int i = 0; i < N_INITIAL_BUFFER; i++) {
        turnIntoBuffer(p[i]);
        p[i].r.y -= W.y;
        particleBuffer.push(i);
    }

    for (int i = 0; i < N_ORIGIN_CELLS; i++) {
        int3 gridDim = INITIAL_ORGANISM_DIMENSIONS; //genomes[iOrigin].gridDim;

        // Define number of in- and outputs
        int inputs = N_INPUTS;              // X, Y, Z, Dist
        int nonCelltypeOutputs = 2;         // Growth prob
        int outputs = N_CELL_TYPES + nonCelltypeOutputs;
        Genome g = Genome(inputs, outputs, gridDim);
        g.mutate(); g.mutate(); g.mutate(); g.mutate(); g.mutate();
        xyz origin = make_xyz_uniform() * int3_to_xyz(W);

        spawnOrganism(fx, pSet, lSet, &nerveMap, origin, &particleBuffer, p, &g, genomes);
    }
    fx->applyParticleArray(pSet);

    for (int step = 0; step < N_STEPS; step++) {
        g.nCells = 0;
        fx->runEach(boundary(), pSet);
        fx->runPair(particlePair(), pSet, pSet, RANGE);
        fx->runLink(transmitNerveSignals(), lSet);
        fx->runEach(nerveActivationFunction(), lSet);
        fx->runLink(updateNerveWeights(), lSet);
        fx->runEach(moveParticle(), pSet);
        fx->runEach(buoyancy(), pSet);
        fx->runEach(handleEnergy(), pSet);
        fx->runEach(integrate(), pSet);
        //parallel_for (int(0), N, [&](int i)
        for (int i = 0; i<N; i++)
        {
            if (p[i].toBuffer) {
                if (particleBuffer.unsafe_size() > N_INITIAL_BUFFER) {
                    turnIntoEnergy(p[i]);                  
                } else {
                    turnIntoBuffer(p[i]);
                    particleBuffer.push(i);
                }
                if (p[i].type == Neuron ||
                    p[i].type == Motor ||
                    p[i].type == Sense)
                {
                    auto range = nerveMap.equal_range(i);
                    for_each(range.first, range.second, [&](intMap::value_type& nerve){
                        fx->removeLink(lSet, pSet, i, pSet, nerve.second);
                        fx->removeLink(lSet, pSet, nerve.second, pSet, i);
                    });
                }
                p[i].toBuffer = false;
                fx->applyParticleArray(pSet);
            }
          if (p[i].particleType == Cell) {
                // Create offspring:
                if (p[i].type == Egg) { //&& p[i].toReproduce) {
                    //printf("Has %.2f, needs %.2f - ", p[i].energy, genomes[i].getMaxCellsReq() * CELL_INITIAL_ENERGY);
                    if (genomes[i].getMaxCellsReq() * CELL_INITIAL_ENERGY <= p[i].energy) {
                        printf("Reproducing!\n");
                        Genome g = genomes[i];
                        g.mutate();
                        spawnOrganism(fx, pSet, lSet, &nerveMap, p[i].r, &particleBuffer, p, &g, genomes);
                        p[i].toReproduce = false;
                        p[i].toBuffer = true;
                        fx->applyParticleArray(pSet);
                    }
                }
            }
        } //);

        if (step % 10 == 0) {
            printf("nCells: %i\t", g.nCells);
            printf("currgenomeIndex: %i\t", currGenomeIndex);
            printf("step %d\n", step);
            fx->outputFrame("output");
        }

//        if (!g.nCells) {
//            break;
//        }
    }
    delete[] genomes;
    delete fx;

    system("shutdown -s -c \"Simulation done, shutting down in two minutes\" -t 120");
}
