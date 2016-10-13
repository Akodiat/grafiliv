//#include "fluidix.h"
#include "C:\Program Files (x86)\Fluidix\include\fluidix.h"
#include "../lib/genome.h"
#include <queue>

#define DT 0.01f // integration time-step
#define CELL_INITIAL_ENERGY 100.0f
#define PELLET_LIFETIME 5.0f

#define W 300
#define N 100000
#define N_ORIGIN_CELLS 2000
#define N_INITIAL_BUFFER 1000
#define N_STEPS 1000000

//Inputs x,y,z,d:
#define N_INPUTS 4

#define RANGE 3.0f
#define MOVE_FACTOR 100

#define REPULSION_FORCE 150
#define ATTRACTION_FORCE 250

#define CELL_MIN_ENERGY 0.5f
#define PELLET_MIN_ENERGY 0.01f
#define DIVISION_ENERGY (CELL_MIN_ENERGY * 8.0f)
#define ENERGY_PARTICLE_ENERGY 10.0f
#define MAX_CELL_DIVISIONS 10

#define CELL_METABOLISM 0.5f
#define CELL_DECAY_RATE 0.1f

#define FLUID_DENSITY 1.0f
#define G -9.81f

#define GREEN   0.5f
#define RED     1.0f
#define YELLOW  0.7f
#define BLUE    0.0f
#define CYAN    0.3f
#define ORANGE  0.8f

#define transmitFloat(a, b, f) {addFloat(a, -f); addFloat(b, f);}

using namespace std;

enum CellType {
    Photo, Digest, Sting, Vascular, Fat, Sense, Motor, Ballast, Egg, 
    N_CELL_TYPES
};
enum ParticleType {
    Cell, Energy, Pellet, Buffert,
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
    float energySharing;
    float maxEnergy;
    int nDivisions;
    float signal;
    int organism;
    bool toGrow;
    bool reproduce;
    Genome genome;
    Particle *origin;
    CellType type;
};

#define turnIntoPellet(p) {         \
    p.particleType = Pellet;        \
    p.energy = PELLET_LIFETIME;     \
    p.density = FLUID_DENSITY * 2;  \
    p.alpha = 0.5f;                 \
    p.organism = NULL;              \
}                                   \

#define turnIntoEnergy(p) {         \
    p.particleType = Energy;        \
    p.r = make_xyz(                 \
        rnd_uniform()*W,            \
        rnd_uniform()*W + W,        \
        rnd_uniform()*W             \
    );                              \
    p.color = 0.7f;                 \
    p.energy = ENERGY_PARTICLE_ENERGY; \
    p.alpha = 0.3f;                 \
    p.radius = 0.5f;                \
    p.density = FLUID_DENSITY * 20; \
}

#define turnIntoBuffert(p) {        \
    p.particleType = Buffert;       \
    p.r = make_xyz(                 \
        rnd_uniform()*W,            \
        rnd_uniform()*W - W,        \
        rnd_uniform()*W             \
    );                              \
    p.density = FLUID_DENSITY * 2;  \
    p.alpha = 0.5f;                 \
    p.color = 0.5f;                 \
    p.radius = 1.0f;                \
    p.organism = NULL;              \
}

FUNC_EACH(init,
    turnIntoEnergy(p);
    //p.r = (make_xyz_uniform() * W) + make_xyz(0, W, 0);
    //p.color = 0.7f;
    //p.alpha = 0.3f;
    //p.radius = 0.5f;
    //p.density = FLUID_DENSITY * 10;
    //p.particleType = Energy;
    //p.toGrow = false;
    //p.growthProb = -1.0f;
    //p.energy = ENERGY_PARTICLE_ENERGY;
    //p.signal = 0.0f;
    //p.organism = -1;
    //p.reproduce = false;
    //p.origin = NULL;
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
            turnIntoBuffert(p)
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
    if (p.particleType == Cell && p.type == Egg) {
        if (rnd_uniform() < p.signal)
            p.reproduce = true;
    }
)

// bouncing hard wall boundary condition
FUNC_EACH(boundary,
    // Check for wierd 1.#R values... NaN?
    if (p.r.x != p.r.x || p.r.y != p.r.y || p.r.z != p.r.z) {
        printf("Weird particle!!! type (%i:%i)\n", p.particleType, p.type);
        turnIntoBuffert(p);
    } if (p.particleType != Buffert) {
        if (p.r.x < 0) { p.v.x = 0.9f * (0 - p.r.x) / DT; p.r.x = 0; }
        if (p.r.x > W) { p.v.x = 0.9f * (W - p.r.x) / DT; p.r.x = W; }
        if (p.r.z < 0) { p.v.z = 0.9f * (0 - p.r.z) / DT; p.r.z = 0; }
        if (p.r.z > W) { p.v.z = 0.9f * (W - p.r.z) / DT; p.r.z = W; }

        if (p.particleType == Energy){
            if (p.r.y < 0) {
                turnIntoBuffert(p);
            }
        }
        else {
            if (p.r.y < 0) { p.v.y = 0.9f * (0 - p.r.y) / DT; p.r.y = 0; }
            if (p.r.y > W) { p.v.y = 0.9f * (W - p.r.y) / DT; p.r.y = W; }
        }
    }
)

FUNC_EACH(growth,
    if (p.particleType == Cell)
        addInteger(g.nCells, 1);

    if (p.particleType == Cell && p.energy > DIVISION_ENERGY && p.nDivisions--) //rnd_uniform() < p.growthProb)
        p.toGrow = true;
    )

#define consumeParticle(a, b) {     \
    addFloat(a.energy, b.energy);   \
    turnIntoBuffert(b);             \
}


//float temp = a.energy;  \
//printf("Particle (type %i-%i) consumes another, energy increases from %.2f to %.2f\n", a.particleType, a.type, temp, a.energy); \

FUNC_PAIR(particlePair,
    if (p1.particleType != Buffert && p2.particleType != Buffert) {
        float ratio = (dr - p1.radius - p2.radius) / (range - p1.radius - p2.radius);
        //float ratio = dr/range;

        float attraction = 0.0f;
        if (p1.particleType == Cell &&
            p2.particleType == Cell)
        {
            //Cells from the same organism
            if (p1.organism == p2.organism) {
                //Attraction between cells of same organism
                attraction = ATTRACTION_FORCE * ratio;

                //Signalling between cells of same organism
                float meanSignal = (p1.signal + p2.signal) / 2;
                p1.signal = p2.signal = meanSignal;

                //Energy transmission
                float p1Surplus = maxf(p1.energy - CELL_MIN_ENERGY, 0);
                float p2Surplus = maxf(p2.energy - CELL_MIN_ENERGY, 0);
                transmitFloat(p1.energy, p2.energy, p1Surplus * p1.energySharing);
                transmitFloat(p2.energy, p1.energy, p2Surplus * p2.energySharing);
            }
            //Cells from different organisms
            else {
                //Damage the other cell if you are sting
                if (p1.type == Sting) p2.energy -= 0.1;
                if (p2.type == Sting) p1.energy -= 0.1;
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

        xyz f = u * (REPULSION_FORCE * (1 - ratio) - attraction);

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
    cell->density = FLUID_DENSITY;
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
        cell->energySharing = 0.8f;
        cell->maxEnergy = 5.0f;
        break;
    case Digest:
        cell->color = RED;
        cell->energySharing = 0.8f;
        cell->maxEnergy = 5.0f;
        break;
    case Fat:
        cell->color = 0.4f;
        cell->energySharing = 0.1f;
        cell->maxEnergy = 20.0f;
        break;
    case Motor:
        cell->color = YELLOW;
        cell->energySharing = 0.0f;
        cell->maxEnergy = 5.0f;
        break;
    case Sense:
        cell->color = BLUE;
        cell->energySharing = 0.5f;
        cell->maxEnergy = 5.0f;
        break;
    case Ballast:
        cell->color = CYAN;
        cell->energySharing = 0.0f;
        cell->maxEnergy = 5.0f;
        break;
    case Egg:
        cell->color = ORANGE;
        cell->energySharing = 0.0f;
        cell->maxEnergy = 10.0f;
        break;
    default:
        cell->energySharing = 0.5f;
        cell->maxEnergy = 5.0f;
    }
    cell->nDivisions = output[N_CELL_TYPES] * MAX_CELL_DIVISIONS;
}

// Initialize new organism, not inheriting anything
void initializeNewOrganism(Particle *cell) {
    // Define number of in- and outputs
    int inputs = N_INPUTS;              // X, Y, Z, Dist
    int nonCelltypeOutputs = 1;         // Growth prob
    int outputs = N_CELL_TYPES + nonCelltypeOutputs;

    cell->genome = Genome(inputs, outputs);
    cell->genome.mutate();

    cell->organism = currGenomeIndex++;
    cell->r = make_xyz(
        rnd_uniform() * W,
        rnd_uniform() * W, // * 0.2f + 0.8f,
        rnd_uniform() * W
        );
    cell->origin = cell;
    cell->reproduce = false;
    setDefaultCellValues(cell);

    vector<float> input(inputs, 0.0f); //Input origin
    vector<float> output = cell->genome.getOutput(input);

    applyPhenotype(output, cell);
}

// Initialize organism, inheriting from parent
void initializeOffspring(Particle *cell) {
    cell->genome.mutate();
    cell->organism = currGenomeIndex++;
    cell->origin = cell;
    cell->reproduce = false;

    // Define number of in- and outputs
    int inputs = N_INPUTS;                      // X, Y, Z, Dist

    vector<float> input(inputs, 0.0f); //Input origin
    vector<float> output = cell->genome.getOutput(input);

    applyPhenotype(output, cell);
}

void growCell(Particle *parent, Particle *child) {
    normal_distribution<float> rndNormal(0.0f, 1.0f);
    
    //Half of parent's energy goes to the child
    parent->energy /= 2;

    //  Copy constructor
    *child = Particle(*parent);

    // Displace particles from each other
    xyz displacement = xyz_norm(
        make_xyz(
        rndNormal(rndGen),
        rndNormal(rndGen),
        rndNormal(rndGen)
        )) * parent->radius;

    parent->r -= displacement;
    child->r += displacement;

    xyz dr = child->r - child->origin->r;

    vector<float> input;
    input.push_back(dr.x);
    input.push_back(dr.y);
    input.push_back(dr.z);
    input.push_back(xyz_len(dr));

    child->genome.mutate();

    vector<float> output = child->genome.getOutput(input);

    //printf("input: ");  for(float i : input)  printf("%.2f ",i); printf("\t");
    //printf("output: "); for(float o : output) printf("%.2f ",o); printf("\n");

    applyPhenotype(output, child);
}

int main() {
    Fluidix<> *fx = new Fluidix<>(&g);
    int setA = fx->createParticleSet(N);

    currGenomeIndex = 0;
    g.nCells = 0;

    fx->runEach(init(), setA);

    Particle *p = fx->getParticleArray(setA);
    for (int i = 0; i < N_ORIGIN_CELLS; i++)
        initializeNewOrganism(&p[i]);
    
    queue<Particle*> particleBuffert;
    for (int i = N_ORIGIN_CELLS; i < N_ORIGIN_CELLS + N_INITIAL_BUFFER; i++) {
        turnIntoBuffert(p[i]);
        p[i].r.y -= W;
        particleBuffert.push(&p[i]);
    }
    fx->applyParticleArray(setA);

    for (int step = 0; step < N_STEPS; step++) {
        g.nCells = 0;
        fx->runEach(boundary(), setA);
        fx->runPair(particlePair(), setA, setA, RANGE);
        fx->runEach(moveParticle(), setA);
        fx->runEach(buoyancy(), setA);
        fx->runEach(integrate(), setA);
        fx->runEach(handleEnergy(), setA);
        fx->runEach(growth(), setA);
        for (int i = 0; i < N; i++) {
            if (p[i].particleType == Buffert) {
                if (particleBuffert.size() > N_INITIAL_BUFFER) {
                    turnIntoEnergy(p[i]);
                    fx->applyParticleArray(setA);
                } else {
                    particleBuffert.push(&p[i]);
                }
            }
            if (p[i].particleType == Cell) {
                if (p[i].toGrow &&
                    p[i].origin != NULL &&
                    p[i].origin->particleType == Cell &&
                    p[i].origin->organism == p[i].organism &&
                    !particleBuffert.empty()
                    ) {
                        Particle *parent = &p[i];
                        Particle *child  = particleBuffert.front();
                        particleBuffert.pop();
                        p[i].toGrow = false;
                        growCell(parent, child);
                        fx->applyParticleArray(setA);
                }
                // Create offspring:
                if (p[i].type == Egg && p[i].reproduce) {
                    initializeOffspring(&p[i]);
                    fx->applyParticleArray(setA);
                }
            }
        }

        if (step % 10 == 0) {
            printf("nCells: %i\t", g.nCells);
            printf("step %d\n", step);
            fx->outputFrame("output");
        }
    }

    delete fx;
}
