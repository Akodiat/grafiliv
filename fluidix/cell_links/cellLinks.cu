//#include "fluidix.h"
#define PARTICLE_BITS 28
#include "C:\Program Files (x86)\Fluidix\include\fluidix.h"
#include "../lib/fluidix_extended.h"
#include "../lib/fluidix_int_all_pairs.h"
#include "../lib/genome.h"
#include <queue>
#include <ppl.h>
#include <concurrent_queue.h>

#define DT 0.01f // integration time-step
#define CELL_INITIAL_ENERGY 10.0f
#define PELLET_LIFETIME 5.0f

#define W make_int3(400, 200, 400)
#define N 100000
#define N_ORIGIN_CELLS 20
#define N_INITIAL_BUFFER 100000
#define N_STEPS 1000000

//Inputs x,y,z,d:
#define N_INPUTS 4

#define RANGE 3.0f
#define MOVE_FACTOR 50

#define REPULSION_FORCE 50

#define CELL_MIN_ENERGY 1.5f
#define PELLET_MIN_ENERGY 0.01f
#define DIVISION_ENERGY (CELL_MIN_ENERGY * 2.5f)
#define ENERGY_PARTICLE_ENERGY 0.5f
#define MAX_CELL_DIVISIONS 100

#define CELL_METABOLISM 0.05f
#define CELL_DECAY_RATE 0.05f

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
    Photo, Digest, Sting, Vascular, Fat, Sense, Motor, Ballast, Egg, 
    N_CELL_TYPES
};
enum ParticleType {
    Cell, Energy, Pellet, Buffer,
    N_PARTICLE_TYPES
};

int currGenomeIndex;

struct Global {
    int nCells;
} g;

struct Particle {
    ParticleType particleType;
    //float growthProb;
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
    bool toGrow;
    bool toBuffer;
    bool toReproduce;
    Particle *origin;
   // Particle *parent;
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

FUNC_EACH(integrate,
    p.v += p.f * DT;
    p.r += p.v * DT;
    p.f = make_xyz(0, 0, 0);
    p.v *= 0.97f;
)

FUNC_EACH(handleEnergy,
    switch (p.particleType) {
    case Cell:
        p.energy -= CELL_METABOLISM * DT;
        if (p.energy < CELL_MIN_ENERGY || p.origin->particleType != Cell)
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
    //p.color = mapf(p.energy, 0.0f, 10.0f, 0.0f, 1.0f);
)

FUNC_EACH(buoyancy,
    if (p.particleType == Cell && p.type == Ballast)
        p.density = clamp(p.density + p.signal, 0.5f, 2.0f);
    float volume = p.radius * p.radius * PI;
    p.f.y += (p.density - FLUID_DENSITY) * G * volume;
)

FUNC_EACH(moveParticle,
    if (p.particleType == Cell && p.type == Motor) {
        xyz f = xyz_norm(p.origin->r - p.r) * MOVE_FACTOR;
        addVector(p.f, f);
    }
)

FUNC_EACH(reproduction,
    if (p.particleType == Cell && p.type == Egg)
        if (p.energy >= p.maxEnergy)
            p.toReproduce = true;
)

// bouncing hard wall boundary condition
FUNC_EACH(boundary,
    // Check for wierd 1.#R values... NaN?
    if (isWeirdParticle(p)) {
        //printf("Weird particle!!! type (%i:%i)\n", p.particleType, p.type);
        turnIntoBuffer(p);
        p.toBuffer = true;
    } if (p.particleType != Buffer) {
        if (p.r.x < 0)   { p.v.x = 0.9f * (0 - p.r.x)   / DT; p.r.x = 0; }
        if (p.r.x > W.x) { p.v.x = 0.9f * (W.x - p.r.x) / DT; p.r.x = W.x; }
        if (p.r.z < 0)   { p.v.z = 0.9f * (0 - p.r.z)   / DT; p.r.z = 0; }
        if (p.r.z > W.z) { p.v.z = 0.9f * (W.z - p.r.z) / DT; p.r.z = W.z; }

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
)

FUNC_EACH(growth,
    if (p.particleType == Cell)
        addInteger(g.nCells, 1);

    if (p.particleType == Cell && p.energy > DIVISION_ENERGY && p.nDivisions--) //rnd_uniform() < p.growthProb)
        p.toGrow = true;
)

FUNC_ALL_PAIRS(linkCells,
    printf("h");
)

inline __host__ __device__ void syncCuda() {
#ifdef __CUDA_ARCH__
    __syncthreads();
#endif
}


/*
#define SPRING_K 1.0f // spring constant
FUNC_EACH(springToParent,
    if (p.particleType == Cell && p.parent != nullptr &&
        p.parent->particleType == Cell &&
        !isWeirdParticlePointer(p.parent)
        ) {
        xyz parentPos = p.parent->r;
        float parentRad = p.parent->radius;
        syncCuda();
        xyz u = (p.r - parentPos);
        float dr = xyz_len(u);
        if (dr > RANGE) {
            turnIntoBuffer(p);
            p.toBuffer = true;
        } else {
            xyz f = u * ((dr - (p.radius + parentRad)) * SPRING_K);
            //printf("parent.r=(%.2f, %.2f, %.2f)\n", p.parent->r.x, p.parent->r.y, p.parent->r.z);
            addVector(p.f, f);
            //addVector(p.parent->f, -f);
        }
    }
)
*/

#define consumeParticle(a, b) {             \
    addFloat(a.energy, b.energy);    \
    b.energy = 0;                           \
    b.toBuffer = true;                      \
}

#define SPRING_K 100.0f // spring constant

FUNC_PAIR(particlePair,
    if (p1.particleType != Buffer && p2.particleType != Buffer) {
        float ratio = (dr - p1.radius - p2.radius) / (range - p1.radius - p2.radius);
        //float ratio = dr/range;

        xyz f = u * (REPULSION_FORCE * (1 - ratio));
        if (p1.particleType == Cell &&
            p2.particleType == Cell)
        {
            //Cells from the same organism
            if (p1.organism == p2.organism) {
                bool neighbours = false;
                for (int n = 0; n < 6; n++)
                if (p1.links[n] == p2_index || p2.links[n] == p1_index) {
                    neighbours = true;
                    break;
                }
                if (neighbours) {
                    //Spring force between neighbours
                    f = -u * ((dr + (p1.radius+p2.radius)) * SPRING_K);

                    //Signalling between cells of same organism
                    float meanSignal = (p1.signal + p2.signal) / 2;
                    p1.signal = p2.signal = meanSignal;

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
                if (p1.type == Sting) turnIntoPellet(p2);
                if (p2.type == Sting) turnIntoPellet(p2);
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

        if (p1.particleType == Cell && p1.type == Sense) p1.signal += 0.1f;
        if (p2.particleType == Cell && p2.type == Sense) p2.signal += 0.1f;

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
    cell->density = FLUID_DENSITY; // * 1.1f;
    cell->particleType = Cell;
}

void applyPhenotype(vector<float> output, Particle *cell) {
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
        cell->energyIn = 0.0f;
        cell->energyOut = 0.2f;
        cell->maxEnergy = 3.0f;
        break;
    case Digest:
        cell->color = RED;
        cell->energyIn = 0.0f;
        cell->energyOut = 0.2f;
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
        cell->energyIn = 0.5f;
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
        cell->maxEnergy = CELL_INITIAL_ENERGY;
        break;
    case Vascular:
        cell->color = 0.2;
        cell->energyIn = 1.0f;
        cell->energyOut = 0.2f;
        cell->maxEnergy = 1.0f;
        break;
    case Sting:
        cell->color = 0.85f;
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
}

#define idxFromCoord(x,y,z) cellBuff.at((x)*yMax*zMax + (y)*zMax + (z))

// Initialize new organism, not inheriting anything
void initializeNewOrganism(xyz origin, vector<int> cellBuff, Particle *p, Genome *genome) {
    int xMax = genome->getGridDim().x;
    int yMax = genome->getGridDim().y;
    int zMax = genome->getGridDim().z;

    int organism = currGenomeIndex++;

    for (int x = 0; x < xMax; x++)
    for (int y = 0; y < yMax; y++)
    for (int z = 0; z < zMax; z++) {
        Particle *cell = &p[idxFromCoord(x, y, z)];
        cell->organism = organism;
        cell->r = origin + make_xyz(x, y, z);
        cell->toReproduce = false;
        setDefaultCellValues(cell);

        vector<float> input;
        input.push_back(x);
        input.push_back(y);
        input.push_back(z);
        input.push_back(xyz_len(make_xyz(x,y,z)));

        vector<float> output = genome->getOutput(input);

        applyPhenotype(output, cell);

        cell->links[0] = x + 1 < xMax ? idxFromCoord(x + 1, y, z) : -1;
        cell->links[1] = y + 1 < yMax ? idxFromCoord(x, y + 1, z) : -1;
        cell->links[2] = z + 1 < zMax ? idxFromCoord(x, y, z + 1) : -1;
        cell->links[3] = x - 1 >= 0 ?   idxFromCoord(x - 1, y, z) : -1;
        cell->links[4] = y - 1 >= 0 ?   idxFromCoord(x, y - 1, z) : -1;
        cell->links[5] = z - 1 >= 0 ?   idxFromCoord(x, y, z - 1) : -1;
    }
}

// Initialize organism, inheriting from parent
void initializeOffspring(Particle *cell, Genome *genome) {
    cell->organism = currGenomeIndex++;
    cell->origin = cell;
//    cell->parent = nullptr;
    cell->toReproduce = false;

    // Define number of in- and outputs
    int inputs = N_INPUTS; // X, Y, Z, Dist

    vector<float> input(inputs, 0.0f); //Input origin
    vector<float> output = genome->getOutput(input);

    applyPhenotype(output, cell);

    printf("New organism! cell type = %i\n", cell->type);
}
/*
void growCell(Particle *parent, Particle *child, Genome *genomeParent, Genome *genomeChild) {
    normal_distribution<float> rndNormal(0.0f, 1.0f);
    
    //Half of parent's energy goes to the child
    parent->energy /= 2;

    //  Copy constructor
    *child = Particle(*parent);
    *genomeChild = Genome(*genomeParent);

    // Displace particles from each other
    xyz displacement = xyz_norm(
        make_xyz(
        rndNormal(rndGen),
        rndNormal(rndGen),
        rndNormal(rndGen)
        )) * parent->radius;

    parent->r -= displacement;
    child->r += displacement;

//    child->parent = parent;

    xyz dr = child->r - child->origin->r;

    vector<float> input;
    input.push_back(dr.x);
    input.push_back(dr.y);
    input.push_back(dr.z);
    input.push_back(xyz_len(dr));

    genomeChild->mutate();

    vector<float> output = genomeChild->getOutput(input);

    //printf("input: ");  for(float i : input)  printf("%.2f ",i); printf("\t");
    //printf("output: "); for(float o : output) printf("%.2f ",o); printf("\n");

    applyPhenotype(output, child);
}
*/

#define printP(chr, p, i) printf("%c\tp[%i].r=(%.2f, %.2f, %.2f)\n", chr, i, p.r.x, p.r.y, p.r.z)

int main() {
    Fluidix<> *fx = new Fluidix<>(&g);
    int setA = fx->createParticleSet(N);

/*    fx->createGlobalArray(&g.link, N);
    fx->getGlobalArray(&g.link);
    for (int i = 0; i < N; i++) {
        g.link[i] = nullptr;
    }
    fx->applyGlobalArray(&g.link);

    for (int i = 0; i < 100; i++){
        printf(g.link[i] == nullptr? "nullptr\n" : " ??? ");
    }
*/
    currGenomeIndex = 0;
    g.nCells = 0;
    Genome *genomes = new Genome[N];
    
    int linkSet = fx->createLinkSet();

    fx->runEach(init(), setA);
    concurrent_queue<int> particleBuffer;
    Particle *p = fx->getParticleArray(setA);

    for (int i = 0; i < N_INITIAL_BUFFER; i++) {
        turnIntoBuffer(p[i]);
        p[i].r.y -= W.y;
        particleBuffer.push(i);
    }

    for (int i = 0; i < N_ORIGIN_CELLS; i++) {
        int3 d = make_int3(10,10,10); //genomes[iOrigin].gridDim;

        // Define number of in- and outputs
        int inputs = N_INPUTS;              // X, Y, Z, Dist
        int nonCelltypeOutputs = 1;         // Growth prob
        int outputs = N_CELL_TYPES + nonCelltypeOutputs;
        Genome g = Genome(inputs, outputs, d);
        g.mutate(); g.mutate(); g.mutate(); g.mutate(); g.mutate();

        int neededParticles = g.getGridDim().x * g.getGridDim().y * g.getGridDim().z;

        if (neededParticles > particleBuffer.unsafe_size()) {
            printf("Not enought particles in buffer\n");
            break;
        }

        vector<int> particles;
        while (neededParticles) {
            int particle;
            if (particleBuffer.try_pop(particle)) {
                particles.push_back(particle);
                neededParticles--;
            }
        }
        xyz origin = make_xyz_uniform() * int3_to_xyz(W);
       
        initializeNewOrganism(origin, particles, p, &g);
    }
    fx->applyParticleArray(setA);

    for (int step = 0; step < N_STEPS; step++) {
        g.nCells = 0;
        fx->runEach(boundary(), setA);
        fx->runPair(particlePair(), setA, setA, RANGE);
  //      fx->runEach(moveParticle(), setA);
        fx->runEach(buoyancy(), setA);
  //      fx->runEach(handleEnergy(), setA);
  //      fx->runEach(growth(), setA);
  //      fx->runEach(reproduction(), setA);
        fx->runEach(integrate(), setA);
        //fx->runAllPairs(linkCells(), setA);
        //parallel_for (int(0), N, [&](int i)
/*        for (int i = 0; i<N; i++)
        {
            if (p[i].toBuffer) {
                if (particleBuffer.unsafe_size() > N_INITIAL_BUFFER) {
                    turnIntoEnergy(p[i]);                  
                } else {
                    turnIntoBuffer(p[i]);
                    particleBuffer.push(i);
                }
                p[i].toBuffer = false;
                fx->applyParticleArray(setA);
            }
            if (p[i].particleType == Cell) {
                // Create offspring:
                if (p[i].type == Egg && p[i].toReproduce) {
                    initializeOffspring(&p[i], &genomes[i]);
                    fx->applyParticleArray(setA);
                }
                // Cell division
                if (p[i].toGrow &&
                    p[i].origin != nullptr &&
                    p[i].origin->particleType == Cell &&
                    p[i].origin->organism == p[i].organism &&
                    !particleBuffer.empty()
                    ) {
                        int parent = i;
                        int child;
                        if (particleBuffer.try_pop(child)) {
                            p[i].toGrow = false;
                            growCell(&p[parent], &p[child], &genomes[parent], &genomes[child]);
                            fx->applyParticleArray(setA);
                        }
                }
            }
        } //);
        */
//        fx->runEach(springToParent(), setA);

        if (step % 1 == 0) {
            printf("nCells: %i\t", g.nCells);
            printf("currgenomeIndex: %i\t", currGenomeIndex);
            printf("step %d\n", step);
            fx->outputFrame("output");
        }

//        if (!g.nCells) {
            /*printf("\nrebooting...\n");
            for (int i = 0; i < N_INITIAL_BUFFER; i++) {
                turnIntoBuffer(p[i]);
                p[i].r.y -= W;
                particleBuffer.push(i);
            }
            initializeNewOrganism(&p[N_INITIAL_BUFFER], &genomes[N_INITIAL_BUFFER]);*/
//            break;
//        }
    }
    delete[] genomes;
    delete fx;

    //system("shutdown -s -c \"Simulation done, shutting down in two minutes\" -t 120");
}
