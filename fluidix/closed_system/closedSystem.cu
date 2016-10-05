#include "fluidix.h"
//#include "C:\Program Files (x86)\Fluidix\include\fluidix.h"
//#include "C:\Users\admin\Documents\Joakim\repo\grafiliv\fluidix\lib\genome.h"
#include "../lib/genome.h"

#define DT 0.01f // integration time-step
#define CELL_LIFETIME 50.0f
#define PELLET_LIFETIME 50.0f

#define BOX_MAX 20
#define BOX_MIN 1

#define W 200
#define N 30000
#define N_ORIGIN_CELLS 500
#define N_STEPS 1000000

//Inputs x,y,z,d:
#define N_INPUTS 4

#define RANGE 3.0f
#define MOVE_FACTOR 200

#define REPULSION_FORCE 150
#define ATTRACTION_FORCE 250

#define FLUID_DENSITY 1.0f
#define G -9.81f

#define GREEN   0.5f
#define RED     1.0f
#define YELLOW  0.7f
#define BLUE    0.0f
#define CYAN    0.3f
#define ORANGE  0.8f

using namespace std;

enum CellType {Photo, Pred, Sense, Move, Ballast, Sex,
	N_CELL_TYPES
};
enum ParticleType {Cell, Energy, Pellet,
	N_PARTICLE_TYPES
};

struct Organism {
	Genome genome;
	int linkSet;
	Organism() {}
	Organism(Genome g): genome(g), linkSet(-1) {}
};

int currGenomeIndex;

struct Global {
} g;

struct Particle {
	ParticleType particleType;
	float growthProb;
	xyz r, v, f;
	float color;
	float radius;
	float alpha;
	float density;
	float lifeTime;
	float signal;
	int organism;
	int toGrow;
	bool reproduce;
	Genome genome;
	Particle *origin;
	CellType type;
};

FUNC_EACH(init,
	p.r = (make_xyz_uniform() * W) + make_xyz(0, W, 0);
	p.color = 0.7f;
	p.alpha = 0.3f;
	p.radius = 0.5f;
	p.density = FLUID_DENSITY * 2;
	p.particleType = Energy;
	p.toGrow = -1;


	p.growthProb = -1.0f;
	p.lifeTime = 0.0f;
	p.signal = 0.0f;
	p.organism = -1;
	p.reproduce = false;
	p.origin = NULL;
)

FUNC_EACH(integrate,
	p.v += p.f * DT;
	p.r += p.v * DT;
	p.f = make_xyz(0, 0, 0);
	p.v *= 0.97f;
)


#define turnIntoPellet(p) {			\
	p.particleType = Pellet;			\
	p.lifeTime = PELLET_LIFETIME;	\
	p.density = FLUID_DENSITY * 2;	\
	p.alpha = 0.5f;						\
	p.organism = NULL;					\
}											\

#define turnIntoEnergy(p) {			\
	p.particleType = Energy;			\
	p.r = make_xyz(						\
		rnd_uniform()*W,				\
		W,									\
		rnd_uniform()*W					\
	);										\
	p.color = 0.7f;						\
	p.alpha = 0.3f;						\
	p.radius = 0.5f;					\
	p.density = FLUID_DENSITY * 2;	\
}

FUNC_EACH(age,
	if (p.particleType != Energy) {
		p.lifeTime -= DT;
		if (p.lifeTime <= 0) {
			if (p.particleType == Cell){
				turnIntoPellet(p);
			} else if (p.particleType == Pellet){
				turnIntoEnergy(p);
			}
		}
		if (p.particleType == Cell && p.origin->particleType != Cell)
			turnIntoPellet(p);
	}
)

FUNC_EACH(buoyancy,
	if(p.particleType == Cell && p.type == Ballast)
		p.density = clamp(p.density + p.signal, 0.5f, 2.0f);
	float volume = p.radius * p.radius * PI;
	p.f.y += (p.density - FLUID_DENSITY) * G * volume;
)

FUNC_EACH(moveParticle,
	if(p.particleType == Cell && p.type == Move) {
		xyz f = xyz_norm(p.origin->r - p.r) * MOVE_FACTOR;
		addVector(p.f, f);
	}
)

FUNC_EACH(reproduce,
	if(p.particleType == Cell && p.type == Sex) {
		if(rnd_uniform() < p.signal)
			p.reproduce = true;
	}
)

// bouncing hard wall boundary condition
FUNC_EACH(boundary,
	if (p.r.x < 0) { p.v.x = 0.9f * (0 - p.r.x) / DT; p.r.x = 0; }
	if (p.r.x > W) { p.v.x = 0.9f * (W - p.r.x) / DT; p.r.x = W; }
	if (p.r.z < 0) { p.v.z = 0.9f * (0 - p.r.z) / DT; p.r.z = 0; }
	if (p.r.z > W) { p.v.z = 0.9f * (W - p.r.z) / DT; p.r.z = W; }

	if(p.particleType == Energy){
		if (p.r.y < 0) {
			p.r.x = rnd_uniform() * W;
			p.r.z = rnd_uniform() * W;
			p.r.y = W;
		}
	} else {
		if (p.r.y < 0) { p.v.y = 0.9f * (0 - p.r.y) / DT; p.r.y = 0; }
		if (p.r.y > W) { p.v.y = 0.9f * (W - p.r.y) / DT; p.r.y = W; }
	}
)

#define harvestParticle(a, b, b_index) { 		\
	if(rnd_uniform() < a.growthProb)				\
		a.toGrow = b_index; 							\
	else {													\
		a.origin->lifeTime	+= CELL_LIFETIME/2;	\
		a.lifeTime 			+= 	CELL_LIFETIME/2;	\
		turnIntoEnergy(b);								\
	}														\
}

FUNC_PAIR(particlePair,
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
			float meanSignal = (p1.signal + p2.signal)/2;
			p1.signal = p2.signal = meanSignal;
		//Cells from differant organisms
		} else {
			//Eat the other cell if you are predatory
			if (p1.type == Pred) harvestParticle(p1, p2, p2_index)
			if (p2.type == Pred) harvestParticle(p2, p1, p1_index)
		}
	//If p1 is a cell
	} else if (p1.particleType == Cell) {
		if ((p1.type == Photo && p2.particleType == Energy) ||
			(p1.type == Pred  && p2.particleType == Pellet)
		) harvestParticle(p1, p2, p2_index)
	}
	//If p2 is a cell
	else if (p2.particleType == Cell) {
		if ((p2.type == Photo && p1.particleType == Energy) ||
			(p2.type == Pred  && p1.particleType == Pellet)
		) harvestParticle(p2, p1, p1_index)
	}

	if(p1.particleType == Cell && p1.type == Sense) p1.signal += 0.1f;
	if(p2.particleType == Cell && p2.type == Sense) p2.signal += 0.1f;

	xyz f = u * (REPULSION_FORCE * (1-ratio) - attraction);

	addVector(p1.f, f);
	addVector(p2.f, -f);

	//p1.color = p1.signal;
	//p2.color = p2.signal;
)

void setDefaultCellValues(Particle *cell) {
	cell->alpha		= 1.0f;
	cell->radius		= 1.0f;
	cell->lifeTime	= CELL_LIFETIME;
	cell->density	= FLUID_DENSITY;
	cell->particleType = Cell;
}

void applyPhenotype(vector<float> output, Particle *cell) {
    float max = output[0]; cell->type = (CellType) 0;
    for(int j=1; j<N_CELL_TYPES; j++) {
        if(output[j] > max) {
            max = output[j];
            cell->type = (CellType) j;
        }
    }
    switch(cell->type) {
        case Photo:     cell->color = GREEN;  break;
        case Pred:      cell->color = RED;    break;
        case Move:      cell->color = YELLOW; break;
        case Sense:     cell->color = BLUE;   break;
        case Ballast:   cell->color = CYAN;   break;
        case Sex:       cell->color = ORANGE; break;
    }
    cell->growthProb = output[N_CELL_TYPES];
}

// Initialize new organism, not inheriting anything
void initializeNewOrganism(Particle *cell) {
	// Define number of in- and outputs
	int inputs = N_INPUTS; 				// X, Y, Z, Dist
	int nonCelltypeOutputs = 1;			// Growth prob
	int outputs = N_CELL_TYPES + nonCelltypeOutputs;

	cell->genome = Genome(inputs, outputs);
	cell->genome.mutate();
	//cell->genome.printMathematica();

	cell->organism	= currGenomeIndex++;
	printf("First batch of species! index: %i\n", cell->organism);
	cell->r			= make_xyz_uniform() * W;
	cell->origin		= cell;
	cell->reproduce	= false;
	setDefaultCellValues(cell);

	vector<float> input(inputs, 0.0f); //Input origin
	vector<float> output = cell->genome.getOutput(input);

   applyPhenotype(output, cell);
}

// Initialize organism, inheriting from parent
void initializeOffspring(Particle *cell) {
	cell->genome.mutate();
	cell->organism = currGenomeIndex++;
	cell->origin		= cell;
	cell->reproduce	= false;
	printf("New species! index: %i\n", cell->organism);

	// Define number of in- and outputs
	int inputs = N_INPUTS; 						// X, Y, Z, Dist

	vector<float> input(inputs, 0.0f); //Input origin
	vector<float> output = cell->genome.getOutput(input);

   applyPhenotype(output, cell);
}

void growCell(Particle *parent, Particle *child) {
	normal_distribution<float> rndNormal(0.0f, 1.0f);

	//	Copy constructor
	*child = Particle(*parent);
	child->lifeTime	= CELL_LIFETIME;

	// Displace particles from each other
	xyz displacement = xyz_norm(
		make_xyz(
			rndNormal(rndGen),
			rndNormal(rndGen),
			rndNormal(rndGen)
		)) * parent->radius;

	parent->r -= displacement;
	child->r  += displacement;

	xyz dr = child->r - child->origin->r;

	vector<float> input;
	input.push_back(dr.x);
	input.push_back(dr.y);
	input.push_back(dr.z);
	input.push_back(xyz_len(dr));

	vector<float> output = child->genome.getOutput(input);

	//printf("input: ");  for(float i : input)  printf("%.2f ",i); printf("\t");
	//printf("output: "); for(float o : output) printf("%.2f ",o); printf("\n");

   applyPhenotype(output, child);
}

int main() {
	Fluidix<> *fx = new Fluidix<>(&g);
	int setA = fx->createParticleSet(N);
	//fx->createGlobalArray(&g.toGrow, N * sizeof(int));
	//for(int i=0; i<N; i++)
	//	g.toGrow[i] = -1;
	currGenomeIndex = 0;

	fx->runEach(init(), setA);

	Particle *p = fx->getParticleArray(setA);
	for(int i=0; i<N_ORIGIN_CELLS; i++)
		initializeNewOrganism(&p[i]);
	fx->applyParticleArray(setA);

	for(int step=0; step < N_STEPS; step++) {
		fx->runPair(particlePair(), setA, setA, RANGE);
		fx->runEach(moveParticle(), setA);
		fx->runEach(buoyancy(), setA);
		fx->runEach(boundary(), setA);
		fx->runEach(integrate(), setA);
		fx->runEach(age(), setA);
		fx->runEach(reproduce(), setA);

		for(int i=0; i<N; i++) {
			if(p[i].particleType == Cell) {
				if(p[i].toGrow != -1) {
					int parent = i;
					int child = p[i].toGrow;
					//printf("toGrow! parent: %i, child: %i\n", parent, child);
					growCell(&p[parent], &p[child]);
					//fx->addLink(
					//	g.organisms[p[i].organism].linkSet,
					//	setA, parent, setA, child
					//);
					p[i].toGrow = -1;
					fx->applyParticleArray(setA);
				}
				// Create offspring:
				else if(p[i].type == Sex && p[i].reproduce) {
					initializeOffspring(&p[i]);
					fx->applyParticleArray(setA);
				}
			}
			// Remove dead organisms:
			//else if(p[i].origin == &p[i])
			//	g.organisms.erase(p[i].organism);
		}

		if (step % 10 == 0) {
			printf("step %d\n", step);
			fx->outputFrame("output");
		}
	}
	//createCells(fx, make_xyz(0, 0, 0), L, setA);

	delete fx;
}
