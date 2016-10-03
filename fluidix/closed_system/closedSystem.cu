#include "fluidix.h"
#include "../lib/genome.h"

#define DT 0.01f // integration time-step
#define CELL_LIFETIME 50.0f
#define PELLET_LIFETIME 50.0f

#define BOX_MAX 20
#define BOX_MIN 1

#define W 200
#define N 40000
#define N_ORIGIN_CELLS 100
#define N_STEPS 100000

#define RANGE 3.0f

#define REPULSION_FORCE 150
#define ATTRACTION_FORCE 250

#define FLUID_DENSITY 1.0f
#define G -9.81f

using namespace std;

enum CellType {Photo, Pred, Sense, Move,
	N_CELL_TYPES
};
enum ParticleType {Cell, Energy, Pellet,
	N_PARTICLE_TYPES
};

struct Organism {
	Genome genome;
	int linkSet;

	Organism(Genome g): genome(g), linkSet(-1) {}
};

struct Global {
	vector<Organism> organisms;
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
	int organism;
	int toGrow;
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
	}
)

FUNC_EACH(buoyancy,
	//if(p.type == Ballast)
	//	p.density += sin(p.signal);
	float volume = p.radius * p.radius * PI;
	p.f.y += (p.density - FLUID_DENSITY) * G * volume;
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
	if(p1.particleType == Cell && 
		p2.particleType == Cell && 
		p1.organism == p2.organism) 
	{
		attraction = ATTRACTION_FORCE * ratio;
	}
	xyz f = u * (REPULSION_FORCE * (1-ratio) - attraction);

	addVector(p1.f, f);
	addVector(p2.f, -f);

	if(	p1.organism != p2.organism && (	// Don't eat yourself.
			p1.particleType == Cell || 		// You have to be a cell to
			p2.particleType == Cell			// be able to eat
	)){
		// Photosynthesis
		if(p1.type == Photo && p2.particleType == Energy) {
			harvestParticle(p1, p2, p2_index)
		} else if(p2.type == Photo && p1.particleType == Energy) {
			harvestParticle(p2, p1, p1_index)
		// Consume other cells or pellets
		} else if(p1.type == Pred && !p2.particleType == Energy) {
			harvestParticle(p1, p2, p2_index)
		} else if(p2.type == Pred && !p1.particleType == Energy) {
			harvestParticle(p2, p1, p1_index)
		}
	}
/*
	if(p1.type == Sense) p1.signal += 0.1f;
	if(p2.type == Sense) p2.signal += 0.1f;

	if(p1.organism == p1.organism){
		float meanSignal = (p1.signal + p2.signal)/2;
		meanSignal = sin(meanSignal);
		p1.signal = p2.signal = meanSignal;
	}
*/
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

// Initialize new organism, not inheriting anything
void initializeNewOrganism(Fluidix<> *fx, Particle *cell, Global &g) {
	// Define number of in- and outputs
	int inputs = 4; 						// X, Y, Z, Dist

	int nonCelltypeOutputs = 1;			// Growth prob
	int outputs = N_CELL_TYPES + nonCelltypeOutputs;

	Organism o(Genome(inputs, outputs));
	o.linkSet = fx->createLinkSet();
	o.genome.mutate();
	o.genome.printMathematica();

	g.organisms.push_back(o);
	cell->organism	= g.organisms.size()-1;
	cell->r			= make_xyz_uniform() * W;
	cell->origin		= cell;
	setDefaultCellValues(cell);

	vector<float> input(inputs, 0.0f); //Input origin
	vector<float> output = o.genome.getOutput(input);

	float max = output[0]; cell->type = (CellType) 0;
	for(int j=1; j<N_CELL_TYPES; j++) {
		if(output[j] > max) {
			max = output[j];
			cell->type = (CellType) j;
		}
	}
	switch(cell->type) {
		case Photo:		cell->color = 0.5f; break; // Green
		case Pred:		cell->color = 1.0f; break; // Red
		case Move:		cell->color = 0.7f; break; // Yellow
		//case Sense:		cell->color = 0.3f; break; // Cyan
		//case Move:		cell->color = 1.0f; break; // Red
		//case Ballast:	cell->color = 0.0f; break; // Blue
	}
	cell->growthProb = output[N_CELL_TYPES];
}

// Initialize organism, inheriting from parent
void initializeOffspring(Fluidix<> *fx, Particle *parent, Particle *child, Global &g) {
	Organism o = Organism(g.organisms[parent->organism]);
	o.linkSet = fx->createLinkSet();
	o.genome.mutate();
	o.genome.printMathematica();

	g.organisms.push_back(o);
	child->organism = g.organisms.size()-1;
}

void growCell(Fluidix<> *fx, Particle *parent, Particle *child) {
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

	Organism o = g.organisms[parent->organism];

	xyz dr = child->r - child->origin->r;

	vector<float> input;
	input.push_back(dr.x);
	input.push_back(dr.y);
	input.push_back(dr.z);
	input.push_back(xyz_len(dr));

	vector<float> output = o.genome.getOutput(input);

	//printf("input: ");  for(float i : input)  printf("%.2f ",i); printf("\t");
	//printf("output: "); for(float o : output) printf("%.2f ",o); printf("\n");

	float max = output[0];
	for(int j=1; j<N_CELL_TYPES; j++) {
		if(output[j] > max) {
			max = output[j];
			child->type = (CellType) j;
		}
	}
	switch(child->type) {
		case Photo:		child->color = 0.5f; break; // Green
		case Pred:		child->color = 1.0f; break; // Red
		case Move:		child->color = 0.7f; break; // Yellow
	}
	child->growthProb = output[N_CELL_TYPES];
}

int main() {
	Fluidix<> *fx = new Fluidix<>(&g);
	int setA = fx->createParticleSet(N);
	//fx->createGlobalArray(&g.toGrow, N * sizeof(int));
	//for(int i=0; i<N; i++)
	//	g.toGrow[i] = -1;

	fx->runEach(init(), setA);
	
	Particle *p = fx->getParticleArray(setA);
	for(int i=0; i<N_ORIGIN_CELLS; i++)
		initializeNewOrganism(fx, &p[i], g);
	fx->applyParticleArray(setA);
	
	for(int step=0; step < N_STEPS; step++) {
		fx->runPair(particlePair(), setA, setA, RANGE);
		fx->runEach(buoyancy(), setA);
		fx->runEach(boundary(), setA);
		fx->runEach(integrate(), setA);
		fx->runEach(age(), setA);
		
		// Grow cells that shold do that
		for(int i=0; i<N; i++) {
			if(p[i].toGrow != -1) {
				int parent = i;
				int child = p[i].toGrow;
				printf("toGrow! parent: %i, child: %i\n", parent, child);
				growCell(fx, &p[parent], &p[child]);
				//fx->addLink(
				//	g.organisms[p[i].organism].linkSet,
				//	setA, parent, setA, child
				//);
				p[i].toGrow = -1;
				fx->applyParticleArray(setA);
			}
		}

		if (step % 10 == 0) {
			printf("step %d\n", step);
			fx->outputFrame("output");
		}
	}
	//createCells(fx, make_xyz(0, 0, 0), L, setA);

	delete fx;
}
