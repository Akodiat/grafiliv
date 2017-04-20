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
int step;

// Turn a particle into Detritus (dead cell) type
// Previous energy is preserved
#define turnIntoDetritus(p) {         \
    p.particleType = Detritus;        \
    p.density = g.fluidDensity * 2; \
    p.organism = -1;                \
}                                   \

// Turn a particle into energy type
#define resetEnergy(p) {               \
    p.particleType = Energy;           \
    p.r.x = rnd_uniform() * g.w.x;     \
    p.r.y += g.w.y;                    \
    p.r.z = rnd_uniform() * g.w.z;     \
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
    p.density = g.fluidDensity;        \
    p.color = 0.5f;                    \
    p.radius = 1.0f;                   \
    p.organism = -1;                   \
}

// Update position r of particles given
// velocity v and force f
FUNC_EACH(integrate,
    p.v += p.f * g.dt;
    p.r += p.v * g.dt;
    p.f = make_xyz(0, 0, 0);
    p.v *= 0.97f;
    )

// Decrease energy through metabolism in
// cells and decay in detritus.
// cell --> detritus --> buffer
FUNC_EACH(handleEnergy,
    p.color = p.energy * 0.2f;
    switch (p.particleType) {
    case Cell:
        //printf("Cell energy: %.2f\n", p.energy);
        p.energy -= p.metabolism * g.dt;
        if (p.energy < g.minCellEnergy)
            turnIntoDetritus(p)
        else if (p.energy > p.maxEnergy)
            p.energy -= (p.energy - p.maxEnergy) * 0.1f;
        break;
    case Detritus:
        p.energy -= g.cellDecayRate * g.dt;
        if (p.energy <= g.minDetritusEnergy)
            p.toBuffer = true;
        break;
    case Energy:
        if (p.energy <= 0)
            p.toBuffer = true;
    }
)

// Particles float depending on their density
FUNC_EACH(buoyancy,
    float volume = sphereVolume(p.radius);
    float weight = p.density * volume;
    float displacedFluidWeight = g.fluidDensity * volume;
    float apparentWeight = weight - displacedFluidWeight;
    p.f.y += apparentWeight * g.gravity;
)

FUNC_EACH(countParticles,
    switch (p.particleType) {
    case Cell:
        addInteger(g.nCells, 1); break;
    case Detritus:
        addInteger(g.nDetritus, 1); break;
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
        printf("Wierd particle at (%.2f, %.2f, %.2f) type=%i\n", p.r.x, p.r.y, p.r.z, p.particleType);
        if (p.particleType == Energy) {
            p.r.y = rnd_uniform() * g.w.y;
            resetEnergy(p);
        }
        else {
            turnIntoBuffer(p);
            p.toBuffer = true;
        }
    } 
    if (p.particleType != Buffer) {
        if (p.r.x < 0)     p.r.x = g.w.x;
        if (p.r.x > g.w.x) p.r.x = 0;
        if (p.r.z < 0)     p.r.z = g.w.z;
        if (p.r.z > g.w.x) p.r.z = 0;

        if (p.particleType == Energy) {
            if (p.r.y < 0) {
                resetEnergy(p);
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

// Let particle a eat of particle b until full
#define getEnergyNeed(a,b) maxf(minf((a.maxEnergy - a.energy), b.energy),0)

// For each particle within a predifined distance
FUNC_PAIR(particlePair,
    if (p1.particleType != Buffer && p2.particleType != Buffer) {
        xyz f = u * maxf(
            (g.repulsiveForce * (1 - dr / (p1.radius + p2.radius))),
            0
        );
        if (p1.particleType == Cell && p2.particleType == Cell)
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
                    //turnIntoDetritus(p2);
                    //transmitFloat(p2.energy, p1.energy, 0.5f);
                    float need = getEnergyNeed(p1, p2);
                    transmitFloat(p2.energy, p1.energy, need);
                }
                if (p2.type == Sting && dr <= (p1.radius + p2.radius)) {
                    //turnIntoDetritus(p2);
                    //transmitFloat(p1.energy, p2.energy, 0.5f);
                    float need = getEnergyNeed(p2, p1);
                    transmitFloat(p1.energy, p2.energy, need);
                }
            }
        }
        //If p1 is a cell
        else if (p1.particleType == Cell && dr <= (p1.radius + p2.radius)) {
            if ((p1.type == Photo && p2.particleType == Energy) ||
                (p1.type == Digest  && p2.particleType == Detritus)
                )
            {
                float need = getEnergyNeed(p1, p2);
                transmitFloat(p2.energy, p1.energy, need);
            }
        }
        //If p2 is a cell
        else if (p2.particleType == Cell && dr <= (p1.radius + p2.radius)) {
            if ((p2.type == Photo && p1.particleType == Energy) ||
                (p2.type == Digest  && p1.particleType == Detritus)
                )
            {
                float need = getEnergyNeed(p2, p1);
                transmitFloat(p1.energy, p2.energy, need);
            }
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
    cell->radius = clamp((radius/2) + 0.5f, 0.5f, 1.0f);
    //if (cell->radius < 0) cell->radius = 0.5f;
    //printf("Radius mapped from %.2f to %.2f\n", radius, cell->radius);
    //float volume = sphereVolume(cell->radius);
    //float mass = 1.0f;

    cell->density = g.fluidDensity * 1.1f; //mass/volume;

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
        cell->energyOut = 0.5f;
        cell->maxEnergy = 10.0f;
        break;
    case Digest:
        cell->energyIn = 0.01f;
        cell->energyOut = 0.5f;
        cell->maxEnergy = 10.0f;
        break;
    case Fat:
        cell->energyIn = 1.0f;
        cell->energyOut = 0.01f;
        cell->maxEnergy = 50.0f;
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
        cell->energyIn = 0.01f;
        cell->energyOut = 0.5f;
        cell->maxEnergy = 10.0f;
        break;
    case Buoyancy:
        cell->energyIn = 1.0f;
        cell->energyOut = 0.0f;
        cell->maxEnergy = 5.0f;
        cell->density = g.fluidDensity * 0.01f;
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
    p[cell].toBuffer = true;
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
    if (nParticlesNeeded > particleBuffer->size()) {
        cerr << "Not enought particles in buffer\n" << endl;
    }
    vector<int> cellBuff;
/*  
    while (nParticlesNeeded) {
        int particle;
        if (particleBuffer->try_pop(particle)) {
            cellBuff.push_back(particle);
            nParticlesNeeded--;
        }
        else
            printf("Failed to retrive from buffer, trying again\n");
    }
*/
    while (nParticlesNeeded--) {
        int particle = particleBuffer->front();
        particleBuffer->pop();
        cellBuff.push_back(particle); 
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
            cell->metabolism += volume * 0.05f;
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

    Organism organism = { genome, nerveSys, cells, -1, 100 };

    //Add organism to organism map
    organisms->emplace(organismID, organism);

    //Output organism to disk
    outputOrganism(&organism, organismID, step);

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
		int3 gridDim = make_int3(1,1,1);

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

    Organism organism = { genome, nerveSys, cells, parent, 100 };

    //Add organism to organism map
    organisms->emplace(organismID, organism);

    //Output organism to disk
    outputOrganism(&organism, organismID, step);

    return organismID;
}

#define printP(chr, p, i) printf("%c\tp[%i].r=(%.2f, %.2f, %.2f)\n", chr, i, p.r.x, p.r.y, p.r.z)

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

    currGenomeIndex = 0;
    OrganismMap organisms;

    // Initialize buffer
    ParticleBuffer particleBuffer;

    Particle *p = fx->getParticleArray(pSet);

    int i = 0;
    int neededEnergy = g.energyParticleCount;
      
    // Initialize energy particles
    while(neededEnergy--) {
        p[i].r.y = rnd_uniform() * g.w.y;
        resetEnergy(p[i]);
        i++;
    }

    // Turn the rest of the particles into buffer
    while (i < g.nParticles) {
        turnIntoBuffer(p[i]);
        particleBuffer.push(i);
        i++;
    }

    //for (int i = 0; i < 1000; i++)
    loadOrg("initOrg.json", &particleBuffer, p, &organisms);

    fx->applyParticleArray(pSet);

    FILE *out = fopen("countCells.csv", "w");
    fprintf(out, "nDetritus,nBuffer,nEnergy,nCells\n");

    step = 0;
    while(step++ < g.nSteps) {
        fx->runEach(boundary(), pSet);
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
                    turnIntoDetritus(p[i]);
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
                    
                    //p[i].f += f;
                    p[i].energy -= xyz_len(f) * g.moveCost;
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

        fx->runEach(buoyancy(), pSet);
        fx->runEach(handleEnergy(), pSet);
        fx->runEach(integrate(), pSet);
        
        //if (step % 100 == 0) { ParticleBuffer empty; swap(particleBuffer, empty); }
        for (int i = 0; i<g.nParticles; i++)
        {
            if (p[i].toBuffer && p[i].particleType != Energy) {
                turnIntoBuffer(p[i]);
                particleBuffer.push(i);
                p[i].toBuffer = false;
                fx->applyParticleArray(pSet);
            }
            //else if (step % 100 == 0 && p[i].particleType == Buffer){
            //    particleBuffer.push(i);
            //}

        }

        g.nDetritus = g.nBuffer = g.nEnergy = g.nCells = 0;
        fx->runEach(countParticles(), pSet);

        //If buffer is getting to small, increase it
        //by adding more particles to the simulation
        if (particleBuffer.size() < g.bufferSize) {
            int currentParticleCount = g.nParticles;
            g.nParticles += g.bufferSize;
            printf("Increasing particle array size from %i to %i\n", currentParticleCount, g.nParticles);
            fx->resizeParticleSet(pSet, g.nParticles);
            p = fx->getParticleArray(pSet);
            for (int i = currentParticleCount; i < g.nParticles; i++) {
                turnIntoBuffer(p[i]);
                particleBuffer.push(i);
            }
            fx->applyParticleArray(pSet);
        }
        else if (particleBuffer.size() > g.nParticles / 2) {
            //printf("Decreasing buffer size from %i", g.nParticles);
            int nBuffersAtEnd = 0;
            //printf("End particle type: %i (buffer is %i)\n", p[g.nParticles - nBuffersAtEnd - 1].particleType, Buffer);
            while (p[g.nParticles - nBuffersAtEnd - 1].particleType == Buffer &&
                particleBuffer.size() - nBuffersAtEnd > g.nParticles / 2
            ){
                nBuffersAtEnd++;
            }
            if (nBuffersAtEnd > 0){
                printf("Decreasing particle array size from %i to %i\n", g.nParticles, g.nParticles - nBuffersAtEnd);
                g.nParticles -= nBuffersAtEnd;
                fx->resizeParticleSet(pSet, g.nParticles);
                p = fx->getParticleArray(pSet);
                fx->applyParticleArray(pSet);

                ParticleBuffer empty;
                swap(particleBuffer, empty);

                for (int i = 0; i<g.nParticles; i++)
                    if (p[i].particleType == Buffer)
                        particleBuffer.push(i);
            }
        }

        if (step % 100 == 0) {
            printf("nOrgs: %i\t", organisms.size());
            printf("currgenomeIndex: %i\t", currGenomeIndex);
            printf("buffer: %i (in queue), %i (actual)\t", particleBuffer.size(), g.nBuffer);
            printf("step %d\n", step);
            outputParticles(p, g.nParticles, step);
            fprintf(out, "%i,%i,%i,%i\n", g.nDetritus, g.nBuffer, g.nEnergy, g.nCells);

        }

        if (step % 10000 == 0) fx->outputFrame("temp");

        if (organisms.size() == 0) {
            printf("All organisms died. End of simulation\n");
            break;
        }
    }
    fclose(out);
    delete fx;
    //system("shutdown -s -c \"Simulation done, shutting down in two minutes\" -t 120");
}
