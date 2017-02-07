//#include "fluidix.h"
#define PARTICLE_BITS 28
#include "C:\Program Files (x86)\Fluidix\include\fluidix.h"
#include "../lib/structures.h"
#include "../lib/genome.h"
#include "../lib/nerveSystem.h"
#include "../lib/io.h"
#include <queue>
#include <ppl.h>

#define transmitFloat(a, b, f) {addFloat(a, -f); addFloat(b, f);}
#define isWeirdParticle(p) (p.r.x != p.r.x || p.r.y != p.r.y || p.r.z != p.r.z)
#define isWeirdParticlePointer(p) (p->r.x != p->r.x || p->r.y != p->r.y || p->r.z != p->r.z)

using namespace std;
using namespace concurrency;


int currGenomeIndex;

#define turnIntoPellet(p) {         \
    p.particleType = Pellet;        \
    p.energy = g.pelletLifetime;    \
    p.density = g.fluidDensity * 2; \
    p.alpha = 0.5f;                 \
    p.organism = -1;                \
}                                   \

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
    p.alpha = 0.3f;                    \
    p.radius = g.energyParticleRadius; \
    p.density = 10.0f;                 \
}

#define turnIntoBuffer(p) {            \
    p.particleType = Buffer;           \
    p.r = make_xyz(                    \
        rnd_uniform() * g.w.x,         \
        rnd_uniform() * g.w.y - g.w.y, \
        rnd_uniform() * g.w.z          \
    );                                 \
    p.density = g.fluidDensity * 2;    \
    p.alpha = 0.1f;                    \
    p.color = 0.5f;                    \
    p.radius = 1.0f;                   \
    p.organism = -1;                   \
}

FUNC_EACH(init,
    turnIntoEnergy(p);
)

FUNC_EACH(integrate,
    p.v += p.f * g.dt;
    p.r += p.v * g.dt;
    p.f = make_xyz(0, 0, 0);
    p.v *= 0.97f;
)

FUNC_EACH(handleEnergy,
    switch (p.particleType) {
    case Cell:
        p.energy -= g.cellMetabolism * g.dt;
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

FUNC_EACH(buoyancy,
    float volume = p.radius * p.radius * PI;
    p.f.y += (p.density - g.fluidDensity) * g.gravity * volume;
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

#define consumeParticle(a, b) {      \
    addFloat(a.energy, b.energy);    \
    b.energy = 0;                    \
    b.toBuffer = true;               \
}

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

FUNC_SURFACE(collideGround,
	if (p.particleType != Energy){
		if (dr > 1) dr = 1;
		p.f += g.groundRepulsiveForce * u * dr;
	}
)

void setDefaultCellValues(Particle *cell) {
    cell->alpha = 1.0f;
    cell->radius = 1.0f;
    cell->energy = g.initialCellEnergy;
    cell->density = g.fluidDensity * 1.10f;
    cell->particleType = Cell;
}

bool applyPhenotype(vector<float> output, Particle *cell) {
    // If cell should not exist, return
    if (output[N_CELL_TYPES] < g.cellExistenceThreshold)
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
        cell->energyIn = 0.01f;
        cell->energyOut = 0.6f;
        cell->maxEnergy = 3.0f;
        break;
    case Digest:
        cell->energyIn = 0.01f;
        cell->energyOut = 0.6f;
        cell->maxEnergy = 3.0f;
        break;
    case Fat:
        cell->energyIn = 1.0f;
        cell->energyOut = 0.01f;
        cell->maxEnergy = 10.0f;
        break;
    case Sense:
        cell->energyIn = 1.0f;
        cell->energyOut = 0.0f;
        cell->maxEnergy = 3.0f;
        break;
    case Egg:
        cell->energyIn = 1.0f;
        cell->energyOut = 0.0f;
        cell->maxEnergy = 1000.0f;
        break;
    case Vascular:
        cell->energyIn = 1.0f;
        cell->energyOut = 0.2f;
        cell->maxEnergy = 1.0f;
        break;
    case Sting:
        cell->energyIn = 1.0f;
        cell->energyOut = 0.0f;
        cell->maxEnergy = 3.0f;
        break;
    default:
        cell->energyIn = 1.0f;
        cell->energyOut = 0.0f;
        cell->maxEnergy = 3.0f;
    }

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

int spawnOrganism(
    xyz origin, ParticleBuffer *particleBuffer,
    Particle *p, Genome genome, NerveSystem nerveSys, OrganismMap *organisms)
{
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
        cell->energy = g.initialCellEnergy;
        setDefaultCellValues(cell);

        vector<float> input;
        input.push_back(x);
        input.push_back(y);
        input.push_back(z);
        input.push_back(xyz_len(make_xyz(x, y, z)));

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

    Organism organism = { genome, nerveSys, addedCells, -1};

    //Add organism to organism map
    organisms->emplace(organismID, organism);

    //Output organism to disk
    outputOrganism(&organism, organismID);

    return organismID;
}

// Initialize new organism
int spawnOrganism(
    xyz origin, ParticleBuffer *particleBuffer,
    Particle *p, int parent, OrganismMap *organisms)
{
	Genome genome;
	
	if(parent == -1){
		int3 gridDim = g.initialOrganismDimensions; //genomes[iOrigin].gridDim;

		// Define number of in- and outputs
		int inputs = g.nGenomeInputs;              // X, Y, Z, Dist
		int nonCelltypeOutputs = 1;         // Cell existence
		int outputs = N_CELL_TYPES + nonCelltypeOutputs;
		genome = Genome(inputs, outputs, gridDim);
		//g.mutate(); g.mutate(); g.mutate(); g.mutate(); g.mutate();
	}
	else {
		genome = Genome(organisms->at(parent).genome);
	}

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
        cell->energy = g.initialCellEnergy;
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

    NerveSystem nervSys;
    if (parent == -1){
        nervSys = NerveSystem(nSensors, 3);
    } else {
        nervSys = NerveSystem(organisms->at(parent).nerveSystem);
        nervSys.updateInputs(nSensors);
	}
	
	nervSys.mutate();
	
    Organism organism = { genome, nervSys, addedCells, parent };
	
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

    float margin = 1.2;
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

int initializeOrganism(ParticleBuffer *particleBuffer, Particle *p, OrganismMap *organisms)
{
    xyz origin = make_xyz_uniform() * int3_to_xyz(g.w);
    origin.y /= 2;
    origin.y += g.w.y / 2;

    return spawnOrganism(origin, particleBuffer, p, -1, organisms);
}

int main() {

    Fluidix<> *fx = new Fluidix<>(&g);

    g = loadConfig("conf.txt");

    int pSet = fx->createParticleSet(g.nParticles);

    //int terrain = generateTerrain(fx);

    currGenomeIndex = 0;
    g.nEggs = 0;
    OrganismMap organisms;

    fx->runEach(init(), pSet);
    ParticleBuffer particleBuffer;
    Particle *p = fx->getParticleArray(pSet);

    int initialBufferSize =
        (g.initialOrganismDimensions.x * 2 + 1) *
        (g.initialOrganismDimensions.y * 2 + 1) *
        (g.initialOrganismDimensions.z * 2 + 1) *
        g.nInitialOrganisms +
        g.bufferSize;

    for (int i = 0; i < initialBufferSize; i++) {
        turnIntoBuffer(p[i]);
        p[i].r.y -= g.w.y;
        particleBuffer.push(i);
    }

    loadOrg("initOrg.json", &particleBuffer, p, &organisms);
    /*
    for (int i = 0; i < g.nInitialOrganisms; i++) {
        initializeOrganism(&particleBuffer, p, &organisms);
    }
    */
    fx->applyParticleArray(pSet);

    int nReboots = 0;
	
    for (int step = 0; step < g.nSteps; step++) {
        g.nEggs = 0;
        fx->runEach(boundary(), pSet);
        //fx->runSurface(collideGround(), terrain, pSet);
        fx->runPair(particlePair(), pSet, pSet, g.interactionRange);


        p = fx->getParticleArray(pSet);
        vector<int> organismsToRemove;
        for (auto& iOrg : organisms) {
            Organism *o = &iOrg.second;
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
            if (nDead > nLiving) {
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
                    p[i].f += f * g.moveFactor;
                    p[i].signal *= 0.5f;
                }
            }
            // Hatch eggs if they have enought energy:
            for (int i : eggs) {
                int maxReqEnergy = 
                    g.initialCellEnergy * 
                    organisms.at(p[i].organism).genome.getMaxCellsReq();

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

        if (step % 10 == 0) {// && step > 1000000
			//|| step % 10000 == 0) {
            printf("nEggs: %i\t", g.nEggs);
            printf("currgenomeIndex: %i\t", currGenomeIndex);
            printf("nReboots: %i\t", nReboots);
            printf("step %d\n", step);
            //if (currGenomeIndex - g.nInitialOrganisms - nReboots > 0)
                outputParticles(p, g.nParticles, step);
        }

        if (!g.nEggs) break;
    }
    delete fx;

    //system("shutdown -s -c \"Simulation done, shutting down in two minutes\" -t 120");
}
