//#include "fluidix.h"
#define PARTICLE_BITS 28
#include "C:\Program Files (x86)\Fluidix\include\fluidix.h"
#include "../lib/genome.h"
#include "../lib/nerveSystem.h"
#include <queue>
#include <ppl.h>
#include <concurrent_queue.h>
#include <unordered_map>

#define DT 0.01f // integration time-step
#define PELLET_LIFETIME 5.0f

#define W make_int3(300, 100, 300)
#define N 300000
#define N_ORIGIN_ORGANISM 1000
#define BUFFER_SIZE 1000
#define N_STEPS 1000000
#define INITIAL_ORGANISM_DIMENSIONS make_int3(1, 1, 1)

//Inputs x,y,z,d:
#define N_INPUTS 4

#define RANGE 10.0f
#define MOVE_FACTOR 100

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


enum CellType {
    Photo, Digest, Sting, Vascular, Fat, Sense, Ballast, Egg,
    N_CELL_TYPES
};
enum ParticleType {
    Cell, Energy, Pellet, Buffer,
    N_PARTICLE_TYPES
};
enum Neigbour { Front, Right, Up, Back, Left, Down };

int currGenomeIndex;

struct Organism {
    Genome genome;
    NerveSystem nervSystem;
    vector<int> cells;
};

typedef std::unordered_map<int, Organism> organismMap;

struct Global {
    int nEggs;
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
    float signal;
    int organism;
    bool toBuffer;
    bool toReproduce;
    int links[6];
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
    p.signal = 0.0f;                \
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
    p.alpha = 0.1f;                 \
    p.color = 0.5f;                 \
    p.radius = 1.0f;                \
    p.organism = -1;                \
}

FUNC_EACH(init,
    turnIntoEnergy(p);
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
/*
FUNC_EACH(moveParticle,
    if (p.particleType == Cell && p.type == Motor) {
        xyz f = make_xyz(0, 0, 0);
        if (p.links[Left]  < 0)  f.x += 1;
        if (p.links[Up]    < 0)  f.y += 1;
        if (p.links[Back]  < 0)  f.z += 1;
        if (p.links[Right] < 0)  f.x -= 1;
        if (p.links[Down]  < 0)  f.y -= 1;
        if (p.links[Front] < 0)  f.z -= 1;
        f = xyz_norm(f) * MOVE_FACTOR * clamp(p.signal, -1.0f, 1.0f);

        addVector(p.f, f);
    }
)
*/

FUNC_EACH(reproduction,
    if (p.particleType == Cell && p.type == Egg)
        if (p.energy >= p.maxEnergy){
            p.toReproduce = true;
            printf("toReproduce!!!\n");
        }
)


// bouncing hard wall boundary condition
FUNC_EACH(boundary,
    if (p.particleType == Cell && p.type == Egg)
        addInteger(g.nEggs, 1);

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
                    if (p1.links[n] == p2_index || p2.links[n] == p1_index) {
                        neighbours = true;
                        break;
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
                    transmitFloat(p2.energy, p1.energy, 0.05f);
                }
                if (p2.type == Sting && dr <= (p1.radius + p2.radius)) {
                    //turnIntoPellet(p2);
                    transmitFloat(p1.energy, p2.energy, 0.05f);
                }
            }
        }
        //If p1 is a cell
        else if (p1.particleType == Cell && dr <= (p1.radius + p2.radius)) {
            if ((p1.type == Photo && p2.particleType == Energy) ||
                (p1.type == Digest  && p2.particleType == Pellet)
                ) consumeParticle(p1, p2)
        }
        //If p2 is a cell
        else if (p2.particleType == Cell && dr <= (p1.radius + p2.radius)) {
            if ((p2.type == Photo && p1.particleType == Energy) ||
                (p2.type == Digest  && p1.particleType == Pellet)
                ) consumeParticle(p2, p1)
        }

        if (p1.particleType == Cell && p1.type == Sense) addFloat(p1.signal, 1.0f);
        if (p2.particleType == Cell && p2.type == Sense) addFloat(p1.signal, 1.0f);

        addVector(p1.f, f);
        addVector(p2.f, -f);

        //p1.color = p1.signal;
        //p2.color = p2.signal;
    }
)

FUNC_SURFACE(collideGround,
if (p.particleType != Energy){
    if (dr > 1) dr = 1;
    p.f += WALL * u * dr;
    /*xyz f = u * 50 * dr;
    p.f += f;
    
    if (p3) f /= 3;
    else if (p2) f /= 2;
    if (p1) addVector(p1->f, -f);
    if (p2) addVector(p2->f, -f);
    if (p3) addVector(p3->f, -f);
    */
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
    case Sense:
        cell->color = BLUE;
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
int spawnOrganism(
    Fluidix<> *fx, int particleSet,
    xyz origin, concurrent_queue<int> *particleBuffer,
    Particle *p, Genome *parentGenome, organismMap *organisms)
{
    Genome genome(*parentGenome);
    genome.mutate();
    
    int nParticlesNeeded = genome.getMaxCellsReq();
    if (nParticlesNeeded > particleBuffer->unsafe_size()) {
        printf("Not enought particles in buffer\n");
        return -1;
    }
    vector<int> cellBuff;
    while (nParticlesNeeded) {
        int particle;
        if (particleBuffer->try_pop(particle)) {
            cellBuff.push_back(particle);
            nParticlesNeeded--;
        }
    }

    int3 br = genome.getBoundingRadius();

    int organismID = currGenomeIndex++;
    vector<int> removedCells;
    vector<int> addedCells;
    
    int nSensors = 0;

    for (int x = -br.x; x <= br.x; x++)
    for (int y = -br.y; y <= br.y; y++)
    for (int z = -br.z; z <= br.z; z++) {
        Particle *cell = &p[iFromCoord(x, y, z)];
        cell->organism = organismID;
        cell->r = origin + make_xyz(x, y, z);
        cell->toReproduce = false;
        cell->energy = CELL_INITIAL_ENERGY;
        setDefaultCellValues(cell);

        vector<float> input;
        input.push_back(x);
        input.push_back(y);
        input.push_back(z);
        input.push_back(xyz_len(make_xyz(x,y,z)));

        vector<float> output = genome.getOutput(input);

        if (applyPhenotype(output, cell)) {
            cell->links[Left] = x + 1 < br.x ? iFromCoord(x + 1, y, z) : -1;
            cell->links[Up] = y + 1 < br.y ? iFromCoord(x, y + 1, z) : -1;
            cell->links[Back] = z + 1 < br.z ? iFromCoord(x, y, z + 1) : -1;
            cell->links[Right] = x - 1 >= 0 ? iFromCoord(x - 1, y, z) : -1;
            cell->links[Down] = y - 1 >= 0 ? iFromCoord(x, y - 1, z) : -1;
            cell->links[Front] = z - 1 >= 0 ? iFromCoord(x, y, z - 1) : -1;

            if (cell->type == Sense)
                nSensors++;
            addedCells.push_back(iFromCoord(x, y, z));
        }
        else
            removedCells.push_back(iFromCoord(x, y, z));
    }
    for (int i : removedCells)
        emptyCellPos(p, i);

    NerveSystem nervSys(nSensors, 3);
    Organism organism = { genome, nervSys, addedCells };
    organisms->emplace(organismID, organism);

    return organismID;
}

#define printP(chr, p, i) printf("%c\tp[%i].r=(%.2f, %.2f, %.2f)\n", chr, i, p.r.x, p.r.y, p.r.z)

int generateTerrain(Fluidix<> *fx){
    exponential_distribution<float> rndUniform(1);

    int terrDimX = 10;
    int terrDimZ = 10;

    int nParticles = (terrDimX*terrDimZ * 2);

    int meshParticles = fx->createParticleSet(nParticles);
    int meshLinks     = fx->createLinkSet();
    Particle *mesh    = fx->getParticleArray(meshParticles);

    float dx = W.x / (terrDimX-1);
    float dz = W.z / (terrDimZ-1);

    for (int x = 0; x < terrDimX; x++)
    for (int z = 0; z < terrDimX; z++){
        int i = x*terrDimZ + z;
        mesh[i].r = make_xyz(
            x*dx,
            rndUniform(rndGen) * 10 + 10,
            z*dz
        );
        mesh[i + nParticles/2].r = make_xyz(
            x*dx,
            0,
            z*dz
        );

        //Link terrain particles together:
        int s = (x - 1)*terrDimZ + z;
        int w = x*terrDimZ + (z - 1);
        int sw = (x - 1)*terrDimZ + (z - 1);

        if ((x - 1) >= 0) fx->addLink(meshLinks, meshParticles, i, meshParticles, s);
        if ((z - 1) >= 0) fx->addLink(meshLinks, meshParticles, i, meshParticles, w);
        if ((x - 1) >= 0 && (z - 1) >= 0) fx->addLink(meshLinks, meshParticles, i, meshParticles, sw);

        if ((x - 1) >= 0) fx->addLink(meshLinks, meshParticles, i + (nParticles / 2), meshParticles, s + (nParticles / 2));
        if ((z - 1) >= 0) fx->addLink(meshLinks, meshParticles, i + (nParticles / 2), meshParticles, w + (nParticles / 2));
        if ((x - 1) >= 0 && (z - 1) >= 0) fx->addLink(meshLinks, meshParticles, i + (nParticles / 2), meshParticles, sw + (nParticles / 2));

        if (x % (terrDimX - 1) == 0 || z % (terrDimZ - 1) == 0){
            fx->addLink(meshLinks, meshParticles, i, meshParticles, i + (nParticles / 2));
            //if (x > 0 && z > 0) fx->addLink(meshLinks, meshParticles, i, meshParticles, s + (nParticles / 2));
        }
        fx->applyParticleArray(meshParticles);
        fx->outputFrame("output");
    }
    return meshLinks;
}

int main() {
    Fluidix<> *fx = new Fluidix<>(&g);
    int pSet = fx->createParticleSet(N);

    int terrain = generateTerrain(fx);

    currGenomeIndex = 0;
    g.nEggs = 0;
    organismMap organisms;

    fx->runEach(init(), pSet);
    concurrent_queue<int> particleBuffer;
    Particle *p = fx->getParticleArray(pSet);

    int initialBufferSize =
        (INITIAL_ORGANISM_DIMENSIONS.x * 2 + 1) *
        (INITIAL_ORGANISM_DIMENSIONS.y * 2 + 1) *
        (INITIAL_ORGANISM_DIMENSIONS.z * 2 + 1) *
        N_ORIGIN_ORGANISM +
        BUFFER_SIZE;

    for (int i = 0; i < initialBufferSize; i++) {
        turnIntoBuffer(p[i]);
        p[i].r.y -= W.y;
        particleBuffer.push(i);
    }

    for (int i = 0; i < N_ORIGIN_ORGANISM; i++) {
        int3 gridDim = INITIAL_ORGANISM_DIMENSIONS; //genomes[iOrigin].gridDim;

        // Define number of in- and outputs
        int inputs = N_INPUTS;              // X, Y, Z, Dist
        int nonCelltypeOutputs = 2;         // Growth prob
        int outputs = N_CELL_TYPES + nonCelltypeOutputs;
        Genome g = Genome(inputs, outputs, gridDim);
        g.mutate(); g.mutate(); g.mutate(); g.mutate(); g.mutate();
        xyz origin = make_xyz_uniform() * int3_to_xyz(W);
        origin.y /= 2;
        origin.y += W.y / 2;

        spawnOrganism(fx, pSet, origin, &particleBuffer, p, &g, &organisms);
    }
    fx->applyParticleArray(pSet);

    for (int step = 0; step < N_STEPS; step++) {
        g.nEggs = 0;
        fx->runEach(boundary(), pSet);
        fx->runSurface(collideGround(), terrain, pSet,RANGE);
        fx->runPair(particlePair(), pSet, pSet, RANGE);


        p = fx->getParticleArray(pSet);
        vector<int> organismsToRemove;
        for (auto& iOrg : organisms) {
            Organism *o = &iOrg.second;
            vector<float> inputs;
            int nLiving = 0;
            int nDead = 0;
            for (int i : o->cells){
                if (p[i].particleType == Cell){
                    if (p[i].type == Sense)
                        inputs.push_back(p[i].signal);
                    nLiving++;
                }
                else
                nDead++;
            }
            if (nDead > nLiving) {
                for (int i : o->cells)
                if (p[i].particleType == Cell)
                    turnIntoPellet(p[i]);
                organismsToRemove.push_back(iOrg.first);
                continue;
            }
            vector<float> output = o->nervSystem.getOutput(inputs);

            xyz f = make_xyz(output[0], output[1], output[2]);
            for (int i : o->cells) {
                p[i].f += f * MOVE_FACTOR;
                p[i].signal *= 0.5f;
            }
        }
        for (int i : organismsToRemove)
            organisms.erase(i);
        fx->applyParticleArray(pSet);

        //fx->runEach(moveParticle(), pSet);
        fx->runEach(buoyancy(), pSet);
        fx->runEach(handleEnergy(), pSet);
        fx->runEach(integrate(), pSet);
        //parallel_for (int(0), N, [&](int i)
        for (int i = 0; i<N; i++)
        {
            if (p[i].toBuffer) {
                if (particleBuffer.unsafe_size() > BUFFER_SIZE) {
                    turnIntoEnergy(p[i]);                  
                } else {
                    turnIntoBuffer(p[i]);
                    particleBuffer.push(i);
                }
                p[i].toBuffer = false;
                fx->applyParticleArray(pSet);
            }
            if (p[i].particleType == Cell) {
                // Create offspring:
                if (p[i].type == Egg) { //&& p[i].toReproduce) {
                    //printf("Has %.2f, needs %.2f - ", p[i].energy, genomes[i].getMaxCellsReq() * CELL_INITIAL_ENERGY);
                    Genome g = organisms.at(p[i].organism).genome;
                    if (g.getMaxCellsReq() * CELL_INITIAL_ENERGY <= p[i].energy) {
                        g.mutate();
                        int orgID = spawnOrganism(fx, pSet, p[i].r, &particleBuffer, p, &g, &organisms);
                        organisms.at(orgID).nervSystem.mutate();
                        p[i].toReproduce = false;
                        p[i].toBuffer = true;
                        fx->applyParticleArray(pSet);
                    }
                }
            }
        }//);

        //fx->applyParticleArray(pSet);

        if (step % 10 == 0) {
            printf("nEggs: %i\t", g.nEggs);
            printf("currgenomeIndex: %i\t", currGenomeIndex);
            printf("step %d\n", step);
            fx->outputFrame("output");
        }

        if (!g.nEggs) {
            break;
        }
    }
    delete fx;

    system("shutdown -s -c \"Simulation done, shutting down in two minutes\" -t 120");
}
