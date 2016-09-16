#include "fluidix.h"

#define W 100 // world size (cube)
#define N 1000 // number of particles
#define STEPS 100000 // number of simulation steps
#define DT 0.01f // integration time-step
#define RANGE 4.0f // fluid interaction cutoff range
#define HARDNESS 20.0f // repulsive force strength
#define FLUID_DENSITY 1.0f

#define REPULSION_FORCE 150
#define ATTRACTION_FORCE 40


#define PHOTOSYNTHESIS_MAX 0.1f
#define CELL_BASIC_METABOLISM 0.3f
#define CELL_BASIC_ENERGY 0.1f
#define CELL_MAX_ENERGY 1.0f
#define CELL_INITIAL_ENERGY 1.0f

#define MUTATION_RATE 0.01f

#define G -9.81f

#define RESURRECT(parent, child) {																							\
	child.energy 					= parent.energy/2;																			\
	parent.energy 				= parent.energy/2;																			\
	child.alpha 					= 1.0f;																						\
	child.remove 					= false;																						\
	child.radius 					= (child.energy + CELL_BASIC_ENERGY)/2;													\
	child.invertSignal 			= rnd_uniform() > MUTATION_RATE ? parent.invertSignal : !parent.invertSignal;	\
	child.density 				= parent.density;																			\
	child.phoSyAbility 			= parent.phoSyAbility 		+ rnd_normal() * MUTATION_RATE;							\
	child.moveAbility 			= parent.moveAbility 		+ rnd_normal() * MUTATION_RATE;							\
	child.senseAbility 			= parent.senseAbility 		+ rnd_normal() * MUTATION_RATE;							\
	child.recieveSigAbility		= parent.recieveSigAbility	+ rnd_normal() * MUTATION_RATE;							\
	child.sendSigAbility 		= parent.sendSigAbility 	+ rnd_normal() * MUTATION_RATE;							\
}

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
	//Node *head;
} g;

/*
#define add_toDivide(val){Node *node = new Node(); node->x = val; node->next = g.head; g.head = node;}

void next_toDivide(){
    Node *n = g.head;
    g.head = g.head->next;
    delete n;
}
bool has_toDivide() {
	return g.head != NULL;
}
*/

// initialize particles
FUNC_EACH(init,

	p.energy = CELL_INITIAL_ENERGY;
	p.r = make_xyz_uniform() * W; 					//Random position within box
	//p.radius = rnd_uniform(); 						//Random size
	p.radius = (p.energy + CELL_BASIC_ENERGY)/2;
	p.density = FLUID_DENSITY; // abs(rnd_normal() / 10 + FLUID_DENSITY); 	//Random density
	p.alpha = 1.0f;
	p.remove = false;
	
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
FUNC_PAIR(pair,
	float ratio = dr / range;
	float interactivity =  (p1.recieveSigAbility * p2.sendSigAbility * p2.recieveSigAbility * p1.sendSigAbility);

	xyz f = u * (REPULSION_FORCE * (1-ratio) * interactivity - ATTRACTION_FORCE * ratio );
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

	if(p1.energy > 0 && p2.energy <= 0){
		RESURRECT(p1, p2);
		//p2.remove = true;
		//addFloat(p1.energy, CELL_BASIC_ENERGY);
	}
	else if (p2.energy > 0 && p1.energy <= 0) {
		RESURRECT(p2, p1);
		//p1.remove = true;
		//addFloat(p2.energy, CELL_BASIC_ENERGY);
	}
	else {
		//Steal energy from other particle proportional to signal alignment with particle direction
		addFloat(p1.energy, xyz_len(xyz_norm(p1.signal) + (-u)) - xyz_len(xyz_norm(p2.signal) + (u)));
		addFloat(p2.energy, xyz_len(xyz_norm(p2.signal) + (u)) - xyz_len(xyz_norm(p1.signal) + (-u)));
		//addFloat(p1.energy, p2.energy/2 - p1.energy/2);
		//addFloat(p2.energy, p1.energy/2 - p2.energy/2);
	} 
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
		p.moveAbility +
		CELL_BASIC_METABOLISM +
		p.senseAbility
	);

	// If dead
	if(p.energy <= 0)
	{
		p.alpha = 0.5f;
		p.density = FLUID_DENSITY;
		p.phoSyAbility = p.moveAbility = p.senseAbility = p.recieveSigAbility = p.sendSigAbility = 0.0f;
	}
	if(p.energy > CELL_MAX_ENERGY){
		p.energy = CELL_MAX_ENERGY;
		//p.energy /= 2;
		//add_toDivide(p_index);
	}

	p.radius = (p.energy + CELL_BASIC_ENERGY)/2;
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

// simulation
int main(int argc, char **argv) {
	Fluidix<> *fx = new Fluidix<>(&g);

	int A = fx->createParticleSet(N);
	fx->runEach(init(), A);

	for (int i = 0; i < STEPS; i++) {
		printf("step %d / %d: ", i, STEPS);
		fx->setTimer();

		// execute interactions
		fx->runPair(pair(), A, A, RANGE);
		fx->runEach(buoyancy(), A);
		fx->runEach(move(), A);
		fx->runEach(integrate(), A);
		fx->runEach(handleEnergy(), A);
		fx->runEach(boundary(), A);
		fx->runEach(deflateSignals(), A);
		fx->runEach(photosynthesis(), A);

		fx->removeParticles(A);

		printf("%.1f ms\n", fx->getTimer());

		// only output to file every 10th step
		if (i % 1 == 0) fx->outputFrame("output");
	}

	delete fx;
}

