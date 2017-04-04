//#include "fluidix.h"
#define PARTICLE_BITS 28
#include "C:\Program Files (x86)\Fluidix\include\fluidix.h"
#include "../lib/structures.h"
#include "../lib/genome.h"
#include "../lib/nerveSystem.h"
#include "../lib/linearAlgebra.h"
#include "../lib/io.h"
#include <queue>
#include <ppl.h>

// Transfer amount f from a to b
#define transmitFloat(a, b, f) {addFloat(a, -f); addFloat(b, f);}

#define sphereVolume(r) (4.0f / 3.0f) * r * r * r * PI; 

// Check if particle position is defined correctly
#define isWeirdParticle(p) (p.r.x != p.r.x || p.r.y != p.r.y || p.r.z != p.r.z)
#define isWeirdParticlePointer(p) (p->r.x != p->r.x || p->r.y != p->r.y || p->r.z != p->r.z)

using namespace std;
using namespace concurrency;

int currGenomeIndex;

// Turn a particle into pellet (dead cell) type
// Previous energy is preserved
#define turnIntoPellet(p) {         \
    p.particleType = Pellet;        \
    p.density = g.fluidDensity * 2; \
    p.organism = -1;                \
}                                   \

// Turn a particle into energy type
#define turnIntoEnergy(p) {            \
    p.particleType = Energy;           \
    p.r = make_xyz(                    \
        rnd_uniform() * g.w.x,         \
        rnd_uniform() * g.w.y + g.w.y, \
        rnd_uniform() * g.w.z          \
    );                                 \
    p.color = 0.7f;                    \
    p.energy = g.energyParticleEnergy; \
    p.signal = 0.0f;                   \
    p.radius = g.energyParticleRadius; \
    p.density = 10.0f;                 \
}

// Turn a pasticle into buffer type
// Place below arena
#define turnIntoBuffer(p) {            \
    p.particleType = Buffer;           \
    p.r = make_xyz(                    \
        rnd_uniform() * g.w.x,         \
        rnd_uniform() * g.w.y - g.w.y, \
        rnd_uniform() * g.w.z          \
    );                                 \
    p.density = g.fluidDensity * 2;    \
    p.color = 0.5f;                    \
    p.radius = 1.0f;                   \
    p.organism = -1;                   \
}

// Initialization function
// Turn all particles into energy particles
FUNC_EACH(init,
    turnIntoEnergy(p);
)

// Update position r of particles given
// velocity v and force f
FUNC_EACH(integrate,
    p.v += p.f * g.dt;
    p.r += p.v * g.dt;
    p.f = make_xyz(0, 0, 0);
    p.v *= 0.97f;
)

// Decrease energy through metabolism in
// cells and decay in pellets.
// cell --> pellet --> buffer
FUNC_EACH(handleEnergy,
    switch (p.particleType) {
    case Cell:
        p.energy -= p.metabolism * g.dt;
        if (p.energy < g.minCellEnergy )
            turnIntoPellet(p)
        else if (p.energy > p.maxEnergy)
            p.energy -= (p.energy - p.maxEnergy) * 0.1f;
        break;
    case Pellet:
        p.energy -= g.cellDecayRate * g.dt;
        if (p.energy <= g.minPelletEnergy)
            p.toBuffer = true;
        break;
    }
)

// Particles float depending on their density
FUNC_EACH(buoyancy,
    float volume = sphereVolume(p.radius);
    p.f.y += (g.fluidDensity - p.density) * volume; //buoyancy
    p.f.y += g.gravity * volume * p.density;        //gravity
)

FUNC_EACH(countParticles,
    switch (p.particleType) {
    case Cell:
        addInteger(g.nCells, 1); break;
    case Pellet:
        addInteger(g.nPellets, 1); break;
    case Buffer:
        addInteger(g.nBuffer, 1); break;
    case Energy:
        addInteger(g.nEnergy, 1); break;
    }
)

// Periodic boundary conditions
FUNC_EACH(boundary,
    // Check for wierd 1.#R values... NaN?
    if (isWeirdParticle(p)) {
        turnIntoBuffer(p);
        p.toBuffer = true;
    } if (p.particleType != Buffer) {
        if (p.r.x < 0)     p.r.x = g.w.x;
        if (p.r.x > g.w.x) p.r.x = 0;
        if (p.r.z < 0)     p.r.z = g.w.z;
        if (p.r.z > g.w.x) p.r.z = 0;

        if (p.particleType == Energy) {
            if (p.r.y < 0) {
                p.toBuffer = true;
            }
        }
        else {
            if (p.r.y < 0) {
				p.v.y = 0.9f * (0 - p.r.y) / g.dt;
				p.r.y = 0;
			}
            if (p.r.y > g.w.y) {
				p.v.y = 0.9f * (g.w.y - p.r.y) / g.dt;
				p.r.y = g.w.y;
			}
        }
    }
    else if (p.r.y < -2 * g.w.y) p.r.y += g.w.y;
)

// Let particle a eat particle b
#define consumeParticle(a, b) {      \
    addFloat(a.energy, b.energy);    \
    b.energy = 0;                    \
    b.toBuffer = true;               \
}

// For each particle within a predifined distance
FUNC_PAIR(particlePair,
    if (p1.particleType != Buffer && p2.particleType != Buffer) {
        xyz f = u * maxf(
            (g.repulsiveForce * (1 - dr / (p1.radius + p2.radius))),
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
                    f = -u * ((dr - (p1.radius + p2.radius)/2) * g.springForce);

                    //Energy transmission
                    float p1Surplus = maxf(p1.energy - g.minCellEnergy, 0);
                    float p2Surplus = maxf(p2.energy - g.minCellEnergy, 0);
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

        if (p1.particleType == Cell && p1.type == Sense) addFloat(p1.signal, 1.0f/dr);
        if (p2.particleType == Cell && p2.type == Sense) addFloat(p1.signal, 1.0f/dr);

        addVector(p1.f, f);
        addVector(p2.f, -f);
    }
)

// Collision with terrain
FUNC_SURFACE(collideGround,
	if (p.particleType != Energy){
		if (dr > 1) dr = 1;
		p.f += g.groundRepulsiveForce * u * dr;
	}
)

// Initialize particle as cell
void setDefaultCellValues(Particle *cell) {
    //cell->radius = 1.0f;
    cell->energy = g.initialCellEnergy;
    //cell->density = g.fluidDensity * 1.10f;
    cell->particleType = Cell;
}

// Given a phenotype network output, apply it to the cell
bool applyPhenotype(vector<float> output, Particle *cell) {
    // If cell should not exist, return
    if (output[N_CELL_TYPES] < g.cellExistenceThreshold)
        return false;
    float radius = output[N_CELL_TYPES + 1];
    cell->radius = clamp(radius, 0.5f, 2.0f);
    float volume = sphereVolume(cell->radius);
    float mass = 1.0f;
    cell->density = mass/volume;

    float max = output[0]; cell->type = (CellType)0;
    for (int j = 1; j<N_CELL_TYPES; j++) {
        if (output[j] > max) {
            max = output[j];
            cell->type = (CellType)j;
        }
    }
    // If no outputs had a positive value, return
    if (max <= 0) {
        return false;
    }
    switch (cell->type) {
    case Photo:
        cell->energyIn = 0.01f;
        cell->energyOut = 0.6f;
        cell->maxEnergy = 5.0f;
        break;
    case Digest:
        cell->energyIn = 0.01f;
        cell->energyOut = 0.6f;
        cell->maxEnergy = 5.0f;
        break;
    case Fat:
        cell->energyIn = 1.0f;
        cell->energyOut = 0.01f;
        cell->maxEnergy = 10.0f;
        break;
    case Sense:
        cell->energyIn = 1.0f;
        cell->energyOut = 0.0f;
        cell->maxEnergy = 5.0f;
        break;
    case Egg:
        cell->energyIn = 1.0f;
        cell->energyOut = 0.0f;
        cell->maxEnergy = 1000.0f;
        break;
    case Vascular:
        cell->energyIn = 1.0f;
        cell->energyOut = 0.2f;
        cell->maxEnergy = 3.0f;
        break;
    case Sting:
        cell->energyIn = 1.0f;
        cell->energyOut = 0.0f;
        cell->maxEnergy = 5.0f;
        break;
    default:
        cell->energyIn = 1.0f;
        cell->energyOut = 0.0f;
        cell->maxEnergy = 5.0f;
    }

    return true;
}

// Remove cell links from and to cell
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

// Helper function to get the 1-dimensional index
// given x,y,z and box size br
int getIdxFromCoord(int x, int y, int z, int3 br)
{
    x += br.x; y += br.y; z += br.z;
    int lY = 2 * br.y + 1;
    int lZ = 2 * br.z + 1;
    return x*lY*lZ + y*lZ + z;
}
#define iFromCoord(x,y,z) cellBuff.at(getIdxFromCoord(x,y,z,br))

// Create cells of an organism given a genome and a nervous system
pair<int, vector<int>> createCellsFromGenotype(
    xyz origin, ParticleBuffer *particleBuffer,
    Particle *p, Genome *genome, NerveSystem *nerveSys, OrganismMap *organisms)
{
    int nParticlesNeeded = genome->getMaxCellsReq();
    if (nParticlesNeeded > particleBuffer->unsafe_size()) {
        cerr << "Not enought particles in buffer\n" << endl;
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
        cell->energy = g.initialCellEnergy;
        cell->metabolism = g.cellMetabolism +
            g.nerveCost * nerveSys->getSize();
        setDefaultCellValues(cell);

        vector<float> input;
        input.push_back(x);
        input.push_back(y);
        input.push_back(z);
        input.push_back(xyz_len(make_xyz(x, y, z)));

        vector<float> output = genome->getOutput(input);

        if (applyPhenotype(output, cell)) {
            cell->links[Left] = x + 1 < br.x ? iFromCoord(x + 1, y, z) : -1;
            cell->links[Up] = y + 1 < br.y ? iFromCoord(x, y + 1, z) : -1;
            cell->links[Back] = z + 1 < br.z ? iFromCoord(x, y, z + 1) : -1;
            cell->links[Right] = x - 1 >= 0 ? iFromCoord(x - 1, y, z) : -1;
            cell->links[Down] = y - 1 >= 0 ? iFromCoord(x, y - 1, z) : -1;
            cell->links[Front] = z - 1 >= 0 ? iFromCoord(x, y, z - 1) : -1;

            if (cell->type == Sense)
                nSensors++;

            float volume = cell->radius * cell->radius * cell->radius * PI * 4 / 3;
            cell->metabolism += volume * 0.05;
            addedCells.push_back(iFromCoord(x, y, z));
        }
        else
            removedCells.push_back(iFromCoord(x, y, z));
    }
    for (int i : removedCells)
        emptyCellPos(p, i);

    nerveSys->updateInputs(nSensors);

    return pair<int, vector<int>>(organismID, addedCells);
}

//Initialize new organism (without parent)
int spawnOrganism(
    xyz origin, ParticleBuffer *particleBuffer,
    Particle *p, Genome genome, NerveSystem nerveSys, OrganismMap *organisms)
{
    pair<int, vector<int>> o = createCellsFromGenotype(
        origin, particleBuffer, p, &genome, &nerveSys, organisms
    );
    int organismID    = o.first;
    vector<int> cells = o.second;

    Organism organism = { genome, nerveSys, cells, -1, 1000 };

    //Add organism to organism map
    organisms->emplace(organismID, organism);

    //Output organism to disk
    outputOrganism(&organism, organismID);

    return organismID;
}

// Initialize new organism from parent
int spawnOrganism(
    xyz origin, ParticleBuffer *particleBuffer,
    Particle *p, int parent, OrganismMap *organisms)
{
	Genome genome;
    NerveSystem nerveSys;

	if(parent == -1) {
		int3 gridDim = g.initialOrganismDimensions;

		// Define number of in- and outputs
		int inputs = g.nGenomeInputs;       // X, Y, Z, Dist
		int nonCelltypeOutputs = 2;         // Cell existence, cell radius
		int outputs = N_CELL_TYPES + nonCelltypeOutputs;
		genome = Genome(inputs, outputs, gridDim);

        int nerveOutputs = 3;
        nerveSys = NerveSystem(nerveOutputs);
	}
	else {
		genome = Genome(organisms->at(parent).genome);
        nerveSys = NerveSystem(organisms->at(parent).nerveSystem);
	}

    genome.mutate();
    nerveSys.mutate();

    pair<int, vector<int>> o = createCellsFromGenotype(
        origin, particleBuffer, p, &genome, &nerveSys, organisms
        );
    int organismID    = o.first;
    vector<int> cells = o.second;

    Organism organism = { genome, nerveSys, cells, parent, 1000 };

    //Add organism to organism map
    organisms->emplace(organismID, organism);

    //Output organism to disk
    outputOrganism(&organism, organismID);

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

    float dx = g.w.x / (terrDimX-1);
    float dz = g.w.z / (terrDimZ-1);

    float margin = 1.2f;
    float shiftX = ((margin - 1)*g.w.x) / 2;
    float shiftZ = ((margin - 1)*g.w.z) / 2;

    for (int x = 0; x < terrDimX; x++)
    for (int z = 0; z < terrDimX; z++){
        int i = x*terrDimZ + z;
        mesh[i].r = make_xyz(
            (x*dx)*margin - shiftX,
            rndUniform(rndGen) * 10 + 10,
            (z*dz)*margin - shiftZ
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

// Initialize a random organism
int initializeOrganism(ParticleBuffer *particleBuffer, Particle *p, OrganismMap *organisms)
{
    xyz origin = make_xyz_uniform() * int3_to_xyz(g.w);
    origin.y /= 2;
    origin.y += g.w.y / 2;

    return spawnOrganism(origin, particleBuffer, p, -1, organisms);
}

Matrix3 getTransform(xyz front, xyz right, xyz up, xyz back, xyz left, xyz down) {
    return Matrix3(
        xyz_norm(right - left),
        xyz_norm(up - down),
        xyz_norm(front - back)
    );
}

int main() {
    // Create Fluidix library object
    Fluidix<> *fx = new Fluidix<>(&g);

    // Load configuration file
    g = loadConfig("conf.txt");

    // Create a particle set with number of particles
    // equal to g.nParticles
    int pSet = fx->createParticleSet(g.nParticles);

    //int terrain = generateTerrain(fx);

    currGenomeIndex = 0;
    OrganismMap organisms;

    // Initialize all particles
    fx->runEach(init(), pSet);

    ParticleBuffer particleBuffer;
    Particle *p = fx->getParticleArray(pSet);

    int initialBufferSize =
        (g.initialOrganismDimensions.x * 2 + 1) *
        (g.initialOrganismDimensions.y * 2 + 1) *
        (g.initialOrganismDimensions.z * 2 + 1) *
        g.nInitialOrganisms +
        g.bufferSize;

    // Turn a large enought number of the particles into buffer
    for (int i = 0; i < initialBufferSize; i++) {
        turnIntoBuffer(p[i]);
        p[i].r.y -= g.w.y;
        particleBuffer.push(i);
    }

    loadOrg("initOrg.json", &particleBuffer, p, &organisms);

    //for (int i = 0; i < g.nInitialOrganisms; i++) {
    //    initializeOrganism(&particleBuffer, p, &organisms);
    //}

    fx->applyParticleArray(pSet);

    FILE *out = fopen("countCells.csv", "w");
    fprintf(out, "nPellets,nBuffer,nEnergy,nCells\n");
    for (int step = 0; step < g.nSteps; step++) {
        fx->runEach(boundary(), pSet);
        //fx->runSurface(collideGround(), terrain, pSet);
        fx->runPair(particlePair(), pSet, pSet, g.interactionRange);

        p = fx->getParticleArray(pSet);
        vector<int> organismsToRemove;
        for (auto& iOrg : organisms) {
            Organism *o = &iOrg.second;

            o->health -= g.dt;

            vector<float> inputs;
            vector<int> eggs;
            int nLiving = 0;
            int nDead = 0;
            for (int i : o->cells){
                if (p[i].particleType == Cell){
                    if (p[i].type == Sense)
                        inputs.push_back(p[i].signal);
                    if (p[i].type == Egg)
                        eggs.push_back(i);
                    nLiving++;
                }
                else
                nDead++;
            }
            if (o->health <= 0 || nDead > nLiving) {
                for (int i : o->cells)
                if (p[i].particleType == Cell)
                    turnIntoPellet(p[i]);
                organismsToRemove.push_back(iOrg.first);
                continue;
            }
            vector<float> output = o->nerveSystem.getOutput(inputs);

            xyz f = make_xyz(output[0], output[1], output[2]);
            for (int i : o->cells) {
                if (p[i].particleType == Cell){
                    xyz front = p[i].links[Front] >= 0 ? p[p[i].links[Front]].r : make_xyz(0, 0, 1);
                    xyz right = p[i].links[Right] >= 0 ? p[p[i].links[Right]].r : make_xyz(1, 0, 0);
                    xyz up = p[i].links[Up] >= 0 ? p[p[i].links[Up]].r : make_xyz(0, 1, 0);
                    xyz back = p[i].links[Back] >= 0 ? p[p[i].links[Back]].r : make_xyz(0, 0, -1);
                    xyz left = p[i].links[Left] >= 0 ? p[p[i].links[Left]].r : make_xyz(-1, 0, 0);
                    xyz down = p[i].links[Down] >= 0 ? p[p[i].links[Down]].r : make_xyz(0, -1, 0);


                    Matrix3 m = getTransform(
                        front - p[i].r,
                        right - p[i].r,
                        up - p[i].r,
                        back - p[i].r,
                        left - p[i].r,
                        down - p[i].r
                    );
                    p[i].f += m.dot(f) * g.moveFactor;
                    //printf("Energy before: %.2f\t", p[i].energy);
                    p[i].energy -= xyz_len(f) * g.moveCost;
                    //printf("Energy after: %.2f\n", p[i].energy);
                    p[i].signal *= 0.5f;
                }
            }
            // Hatch eggs if they have enought energy:
            for (int i : eggs) {
                int maxReqEnergy =
                    g.initialCellEnergy *
                    o->genome.getMaxCellsReq() +
                    o->genome.getSize() * g.genomeCost;

                if (p[i].energy >= maxReqEnergy + g.initialCellEnergy) {
                    spawnOrganism(
                        p[i].r, &particleBuffer,
                        p, p[i].organism, &organisms
                    );
                    p[i].energy -= maxReqEnergy;
                    fx->applyParticleArray(pSet);
                }
            }
        }
        for (int i : organismsToRemove)
            organisms.erase(i);
        fx->applyParticleArray(pSet);

        //fx->runEach(moveParticle(), pSet);
        fx->runEach(buoyancy(), pSet);
        fx->runEach(handleEnergy(), pSet);
        fx->runEach(integrate(), pSet);
        parallel_for (int(0), g.nParticles, [&](int i)
        //for (int i = 0; i<g.nParticles; i++)
        {
            if (p[i].toBuffer) {
                if (particleBuffer.unsafe_size() > g.bufferSize) {
                    turnIntoEnergy(p[i]);
                } else {
                    turnIntoBuffer(p[i]);
                    particleBuffer.push(i);
                }
                p[i].toBuffer = false;
                fx->applyParticleArray(pSet);
            }
        });

        //fx->applyParticleArray(pSet);

        g.nPellets = g.nBuffer = g.nEnergy = g.nCells = 0;
        fx->runEach(countParticles(), pSet);
        fprintf(out, "%i,%i,%i,%i\n", g.nPellets, g.nBuffer, g.nEnergy, g.nCells);

        if (step % 10 == 0) {// && step > 1000000
			//|| step % 10000 == 0) {
            printf("nOrgs: %i\t", organisms.size());
            printf("currgenomeIndex: %i\t", currGenomeIndex);
            printf("buffer: %i (%i)\t", particleBuffer.unsafe_size(), g.nBuffer);
            printf("step %d\n", step);
            outputParticles(p, g.nParticles, step);
        }

        if (organisms.size() == 0) {
            printf("All organisms died. End of simulation\n");
            break;
        }
    }
    fclose(out);
    delete fx;

    //system("shutdown -s -c \"Simulation done, shutting down in two minutes\" -t 120");
}
