#include "fluidix.h"
#include "../lib/genome.h"

#define DT 0.1f // integration time-step

#define BOX_MAX 20
#define BOX_MIN 1

#define W 100

using namespace std;

enum CellType {Blank, Sense, Move, Ballast};

struct Global {
	float grid_pack;
	int3 grid_num;
	xyz grid_origin;
	vector<Genome> organisms;
} g;

struct Particle {
	bool isCell;
	xyz r, v, f;
	float color;
	float radius;
	float alpha;
	int organism;
	Particle *origin;
	CellType type;
};

FUNC_EACH(init,
	p.r = (make_xyz_uniform() * W) + make_xyz(0 , W, 0);
	p.color = 0.7f;
	p.alpha = 0.5f;
	p.radius = 0.5f;
	p.isCell = false;
)

FUNC_EACH(integrate,
	p.v += p.f * DT;
	p.r += p.v * DT;
	p.f = make_xyz(0, 0, 0);
	p.v *= 0.97f;
)

FUNC_EACH(gravity,
	p.f.y -= 9.81f;
)

// bouncing hard wall boundary condition
FUNC_EACH(boundary,
	if (p.r.x < 0) { p.v.x = 0.9f * (0 - p.r.x) / DT; p.r.x = 0; }
	if (p.r.x > W) { p.v.x = 0.9f * (W - p.r.x) / DT; p.r.x = W; }
	if (p.r.z < 0) { p.v.z = 0.9f * (0 - p.r.z) / DT; p.r.z = 0; }
	if (p.r.z > W) { p.v.z = 0.9f * (W - p.r.z) / DT; p.r.z = W; }

	if(p.isCell){
		if (p.r.y < 0) { p.v.y = 0.9f * (0 - p.r.y) / DT; p.r.y = 0; }
		if (p.r.y > W) { p.v.y = 0.9f * (W - p.r.y) / DT; p.r.y = W; }
	} else {
		if (p.r.y < 0) p.r.y = W;
	}
)

// Initialize new organism, not inheriting anything
void initializeNewOrganism(Particle *cell, Global &g){
	int inputs = 4;
	int outputs = 4;

	Genome genome(inputs, outputs);
	genome.printMathematica();

	g.organisms.push_back(genome);
	cell->organism	= g.organisms.size()-1;
	cell->alpha		= 1.0f;
	cell->radius		= 1.0f;
	cell->r			= make_xyz_uniform() * W;
	cell->isCell		= true;

	vector<float> input(inputs, 0.0f);
	vector<float> output = genome.getOutput(input);

	float max = output[0];
	for(int j=1; j<output.size(); j++) {
		if(output[j] > max) {
			max = output[j];
			cell->type = (CellType) j;
		}
	}
	switch(cell->type) {
		case Blank:		cell->color = 0.6f; break;
		case Sense:		cell->color = 0.3f; break;
		case Move:		cell->color = 1.0f; break;
		case Ballast:	cell->color = 0.0f; break;
	}
}

// Initialize organism, inheriting from parent
void initializeOffspring(Particle *parent, Particle *child, Global &g){
	Genome genome = g.organisms[parent->organism];
	genome.mutate();
	genome.printMathematica();

	g.organisms.push_back(genome);
	child->organism = g.organisms.size()-1;
}

void growCell(Particle *parent, Particle *child) {
	//	Copy constructor?
	*child = Particle(*parent);

	// Displace particles from each other
	xyz displacement = xyz_norm(make_xyz_uniform()) * parent->radius;		
	parent->r -= displacement;
	child->r  += displacement;

	Genome genome = g.organisms[parent->organism];

	xyz dr = parent->r - parent->origin->r;

	vector<float> input;
	input.push_back(dr.x);
	input.push_back(dr.y);
	input.push_back(dr.z);
	input.push_back(xyz_len(dr));

	vector<float> output = genome.getOutput(input);

	//printf("input: ");  for(float i : input)  printf("%.2f ",i); printf("\t");
	//printf("output: "); for(float o : output) printf("%.2f ",o); printf("\n");

	float max = output[0];
	for(int j=1; j<output.size(); j++) {
		if(output[j] > max) {
			max = output[j];
			child->type = (CellType) j;
		}
	}
	switch(child->type) {
		case Blank:		child->color = 0.6f; break;
		case Sense:		child->color = 0.3f; break;
		case Move:		child->color = 1.0f; break;
		case Ballast:	child->color = 0.0f; break;
	}
}

int main() {
	Fluidix<> *fx = new Fluidix<>(&g);
	int setA = fx->createParticleSet(10000);
	fx->runEach(init(), setA);
	
	Particle *p = fx->getParticleArray(setA);
	for(int i=0; i<10; i++) {
		initializeNewOrganism(&p[i], g);
		fx->applyParticleArray(setA);
	}
	
	for(int step=0; step < 10000; step++) {
		fx->runEach(integrate(), setA);
		fx->runEach(gravity(), setA);
		fx->runEach(boundary(), setA);
		fx->outputFrame("output");
	}
	//createCells(fx, make_xyz(0, 0, 0), L, setA);

	delete fx;
}
