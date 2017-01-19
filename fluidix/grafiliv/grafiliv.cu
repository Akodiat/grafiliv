#include "fluidix.h"
#include "../lib/genome.h"
#include <map>

#define W 500 // world size (cube)
#define N 500 // number of particles
#define STEPS 100000 // number of simulation steps
#define DT 0.01f // integration time-step
#define RANGE 4.0f // fluid interaction cutoff range
#define HARDNESS 20.0f // repulsive force strength
#define FLUID_DENSITY 1.0f

#define REPULSION_FORCE 150
#define ATTRACTION_FORCE 40


#define PHOTOSYNTHESIS_MAX 1.0f
#define CELL_BASIC_METABOLISM 0.3f
#define CELL_BASIC_ENERGY 2.0f
#define CELL_MAX_ENERGY 1.0f
#define CELL_INITIAL_ENERGY 4000.0f

#define MUTATION_RATE 0.01f

#define G -9.81f

struct Particle {
	xyz r, v, f; 				// position, velocity, force
	xyz signal;
	xyz sense;
	float density;
	float radius;
	float energy;
	float color; 				// display color
	float alpha;
	bool remove;
	bool divide;
	int linkSet;
	int origin;

	//Genome:
	float phoSyAbility; 		// Ability to photosynthesize
	float moveAbility;			// Ability to move
	float senseAbility;			// Ability to sense surroundings
	bool  invertSignal;			// Inverts signal before acting or transmitting
	float recieveSigAbility;  	// Ability to recieve signals	
	float sendSigAbility;		// Ability to send signals
};

struct Node {
    int x;
    Node *next;
};

struct Global {
} g;

// initialize particles
FUNC_EACH(init,

	p.energy = CELL_INITIAL_ENERGY;
	p.r = make_xyz_uniform() * W; 					//Random position within box
	p.radius = 1.0f; //(p.energy + CELL_BASIC_ENERGY)/2;
	p.density = FLUID_DENSITY; // abs(rnd_normal() / 10 + FLUID_DENSITY); 	//Random density
	p.alpha = 1.0f;
	p.remove = false;	
	p.divide = false;
	
	p.invertSignal 		= rnd_uniform() < 0.5f;
	p.phoSyAbility 		= 0.0f; // abs(rnd_normal());
	p.moveAbility 		= 0.0f; // abs(rnd_normal());
	p.senseAbility 		= 0.0f; // abs(rnd_normal());
	p.recieveSigAbility = 0.0f; // abs(rnd_normal());
	p.sendSigAbility 	= 0.0f; // abs(rnd_normal());

	//p.color = p.moveAbility;
	
	p.sense  = make_xyz(0, 0, 0);
	p.signal = make_xyz_uniform()*5;	//Random initial signal
)

// linear repulsion + attraction at distance
FUNC_PAIR(particlePair,
	float ratio = dr / range;
	float attraction = 0.0f;
	if(p1.linkSet == p2.linkSet)
		attraction = ATTRACTION_FORCE * ratio;
	xyz f = u * (REPULSION_FORCE * (1-ratio) - attraction);
	addVector(p1.f, f);
	addVector(p2.f, -f);

	// Update sense values
	addVector(p1.sense, -u);	// These should be normalized...
	addVector(p2.sense, u);

	// Update signal values with sense
	addVector(p1.signal, p1.sense * p1.senseAbility);
	addVector(p2.signal, p2.sense * p2.senseAbility);

	// Transmit signal
	addVector(p1.signal, p2.signal * p1.recieveSigAbility * p2.sendSigAbility * (p1.invertSignal ? -1 : 1));
	addVector(p2.signal, p1.signal * p2.recieveSigAbility * p1.sendSigAbility * (p2.invertSignal ? -1 : 1));

	// p1 consumes p2
	if(p1.energy > 0 && p2.energy <= 0){
		p2.remove = true;
		addFloat(p1.energy, CELL_BASIC_ENERGY);
	}
	// p2 consumes p1
	else if (p2.energy > 0 && p1.energy <= 0) {
		p1.remove = true;
		addFloat(p2.energy, CELL_BASIC_ENERGY);
	}
/*
	else {
		//Steal energy from other particle proportional to signal alignment with particle direction
		addFloat(p1.energy, xyz_len(xyz_norm(p1.signal) + (-u)) - xyz_len(xyz_norm(p2.signal) + (u)));
		addFloat(p2.energy, xyz_len(xyz_norm(p2.signal) + (u)) - xyz_len(xyz_norm(p1.signal) + (-u)));
		//addFloat(p1.energy, p2.energy/2 - p1.energy/2);
		//addFloat(p2.energy, p1.energy/2 - p2.energy/2);
	}
*/

)

// buoyancy 
FUNC_EACH(buoyancy,
	float volume = p.radius * p.radius * PI;
	p.f.y += (p.density - FLUID_DENSITY) * G * volume;
)

// photosynthesis 
FUNC_EACH(photosynthesis,
	p.energy += PHOTOSYNTHESIS_MAX * p.phoSyAbility;
)

// move 
FUNC_EACH(move,
	xyz f = p.signal * p.moveAbility * 0.1f;
	addVector(p.f, f);
)

// handle energy usage
FUNC_EACH(handleEnergy,
	p.energy -= (
		CELL_BASIC_METABOLISM +
		p.moveAbility +
		p.senseAbility
	);

	// If dead
	if(p.energy <= 0)
	{
		p.alpha = 0.5f;
		p.density = FLUID_DENSITY;
		p.phoSyAbility = p.moveAbility = p.senseAbility = p.recieveSigAbility = p.sendSigAbility = 0.0f;
	}
	// If there is enought energy to divide
	if(p.energy > CELL_MAX_ENERGY){
		p.divide = true;
	}

	//p.radius = (p.energy + CELL_BASIC_ENERGY)/2;
)

// Euler integration
FUNC_EACH(integrate,
	//p.color = xyz_len(p.f) / 50.0f;

	p.v += p.f * DT;
	p.r += p.v * DT;
	p.f = make_xyz(0, 0, 0);
)

// bouncing hard wall boundary condition
FUNC_EACH(boundary,
	if (p.r.x < 0) { p.v.x = 0.9f * (0 - p.r.x) / DT; p.r.x = 0; }
	if (p.r.x > W) { p.v.x = 0.9f * (W - p.r.x) / DT; p.r.x = W; }
	if (p.r.y < 0) { p.v.y = 0.9f * (0 - p.r.y) / DT; p.r.y = 0; }
	if (p.r.y > W) { p.v.y = 0.9f * (W - p.r.y) / DT; p.r.y = W; }
	if (p.r.z < 0) { p.v.z = 0.9f * (0 - p.r.z) / DT; p.r.z = 0; }
	if (p.r.z > W) { p.v.z = 0.9f * (W - p.r.z) / DT; p.r.z = W; }
)

// deflate old signals so that they won't overflow 
FUNC_EACH(deflateSignals,
	//p.color = xyz_len(p.signal)/100;
	p.color = p.senseAbility;
	p.signal = xyz_norm(p.signal);
	p.sense = make_xyz(0, 0, 0);
)

int main(int argc, char **argv) {
	Fluidix<> *fx = new Fluidix<>(&g);

	map<int, Genome> genomes;

	int A = fx->createParticleSet(10);
	fx->runEach(init(), A);
	
	Particle *pArray = fx->getParticleArray(A);
	for (int step = 0; step < 5000; step++) {
		
		Particle *p = fx->getParticleArray(A);
		int n = fx->getParticleCount(A);
		for (int i = 0; i < n; i++) {
			if (p[i].remove == true)
			{
				//printf("Removing particle (%d)\n", i);
/*				memcpy(&p[i], &p[n-1], sizeof(Particle)); 	// replace current with last
				fx->applyParticleArray(A); 					// apply changes before any operation
				fx->resizeParticleSet(A, --n); 				// delete the last particle and decrease n
				Particle *p = fx->getParticleArray(A);
*/			}
			else if (p[i].divide == true)
			{
				//printf("Dividing particle (%d)\n", i);
				
				// If this is the origin cell
				if(p[i].linkSet == -1) {
					p[i].linkSet = fx->createLinkSet();

					int inputs = 4;
					int outputs = 3;

					Genome g(inputs, outputs);

					g.mutate();
					g.mutate();
					g.mutate();
					g.mutate();

					g.printMathematica();

					genomes.insert({p[i].linkSet,	g});

					vector<float> input(inputs, 0.0f);
					vector<float> output = g.getOutput(input);
	
					p[i].color = output[0];
					p[i].density = output[1];
					p[i].remove = output[2]<=0;
				}

				p[i].divide = false;
				p[i].energy -= 0.1f;
				p[i].energy /= 2;
				fx->applyParticleArray(A); 
				fx->resizeParticleSet(A, ++n); 				// add a new particle and increase n
				Particle *p = fx->getParticleArray(A);

				memcpy(&p[n-1], &p[i], sizeof(Particle)); 	// copy current to last

				// Displace particles from each other
				xyz dr = xyz_norm(make_xyz_uniform()) * p[i].radius;		
				p[i].r 	-= dr;
				p[n-1].r 	+= dr;

				//Create link between new particle and origin
				//fx->addLink(p[n-1].linkSet, A, p[n-1].origin, A, n-1);
				printf("parent link set: %i\tchild link set: %i\t",p[i].linkSet, p[n-1].linkSet);
				fx->addLink(p[n-1].linkSet, A, i, A, n-1);

				Genome genome = genomes.at(p[n-1].linkSet);

				dr = p[n-1].r - p[p[n-1].origin].r;

				vector<float> input;
				input.push_back(dr.x);
				input.push_back(dr.y);
				input.push_back(dr.z);
				input.push_back(xyz_len(dr));

				vector<float> output = genome.getOutput(input);

				printf("input: ");  for(float i : input)  printf("%.2f ",i); printf("\t");
				printf("output: "); for(float o : output) printf("%.2f ",o); printf("\n");

				p[n-1].color = output[0];
				p[n-1].density = output[1];
				p[n-1].remove = output[2]<=0;
	
				fx->applyParticleArray(A); 					// apply changes before any operation
			}
		}
		//printf("number of particles: %d\n", n);

		fx->runPair(particlePair(), A, A, RANGE);
		fx->runEach(integrate(), A);
		fx->runEach(buoyancy(), A);
		//fx->runEach(photosynthesis(), A);	
		//fx->runEach(friction(), A);
		fx->runEach(boundary(), A);

		if (step % 1 == 0) {
			//printf("step %d\n", step);
			fx->outputFrame("output");
		}
	}

	delete fx;
}

