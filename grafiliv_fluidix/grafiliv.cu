#include "fluidix.h"

#define W 100 // world size (cube)
#define N 1000 // number of particles
#define STEPS 100000 // number of simulation steps
#define DT 0.01f // integration time-step
#define RANGE 2.0f // fluid interaction cutoff range
#define HARDNESS 20.0f // repulsive force strength
#define FLUID_DENSITY 1.0f

#define REPULSION_FORCE 150
#define ATTRACTION_FORCE 40


#define PHOTOSYNTHESIS_MAX 10.0f

#define G -9.81f

struct Global {
} g;

struct Particle {
	xyz r, v, f; 				// position, velocity, force
	xyz signal;
	xyz sense;
	float density;
	float radius;
	float energy;
	float color; 				// display color
	float alpha;

	//Genome:

	float phoSyAbility; 		// Ability to photosynthesize
	float moveAbility;			// Ability to move
	float senseAbility;			// Ability to sense surroundings

	bool  invertSignal;			// Inverts signal before acting or transmitting
	float recieveSigAbility;  	// Ability to recieve signals	
	float sendSigAbility;		// Ability to send signals
};

// initialize particles
FUNC_EACH(init,

	p.energy = 10.0f;
	p.r = make_xyz_uniform() * W; 					//Random position within box
	p.radius = abs(rnd_normal()); 					//Random size
	p.density = abs(rnd_normal()+FLUID_DENSITY); 	//Random density
	p.alpha = 1.0f;
	
	p.invertSignal 		= rnd_uniform() < 0.5f;
	p.phoSyAbility 		= abs(rnd_normal());
	p.moveAbility 		= abs(rnd_normal());
	p.senseAbility 		= abs(rnd_normal());
	p.recieveSigAbility = abs(rnd_normal());
	p.sendSigAbility 	= abs(rnd_normal());

	//p.color = p.moveAbility;
	
	
	p.sense  = make_xyz(0, 0, 0);
	p.signal = make_xyz_uniform()*5;	//Random initial signal
)

// linear repulsion + attraction at distance
FUNC_PAIR(pair,
	float ratio = dr / range;

	xyz f = u * (REPULSION_FORCE * (1-ratio) - ATTRACTION_FORCE * ratio);
	addVector(p1.f, f);
	addVector(p2.f, -f);

	// Update sense values
	addVector(p1.sense, 2*(-u));	// These should be normalized...
	addVector(p2.sense, 2*u);

	// Update signal values with sense
	addVector(p1.signal, p1.sense * p1.senseAbility);
	addVector(p2.signal, p2.sense * p2.senseAbility);

	// Transmit signal
	addVector(p1.signal, p2.signal * p1.recieveSigAbility * p2.sendSigAbility * (p1.invertSignal ? -1 : 1));
	addVector(p2.signal, p1.signal * p2.recieveSigAbility * p1.sendSigAbility * (p2.invertSignal ? -1 : 1));
	
	p1.color = xyz_len(p1.signal)/100;
	p2.color = xyz_len(p2.signal)/100;
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
	xyz f = p.signal * p.moveAbility * 0.001f;
	addVector(p.f, f);
)

// handle energy usage
FUNC_EACH(handleEnergy,
	p.energy -= (
		p.moveAbility * 5 + 
		p.senseAbility * 3
	);
	//p.color = p.energy;

	// If dead
	if(p.energy <= 0)
	{
		p.alpha = 0.5f;
		p.density = FLUID_DENSITY;
		p.phoSyAbility = p.moveAbility = p.senseAbility = p.recieveSigAbility = p.sendSigAbility = 0.0f;
	}
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
	p.signal = xyz_norm(p.signal);
	p.sense = xyz_norm(p.signal);
	
	//p.color = p.signal;
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

		printf("%.1f ms\n", fx->getTimer());

		// only output to file every 10th step
		if (i % 10 == 0) fx->outputFrame("sample");
	}

	delete fx;
}

