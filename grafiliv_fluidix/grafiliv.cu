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

#define G -9.81f

struct Global {
} g;

struct Particle {
	xyz r, v, f; 		// position, velocity, force
	float density;
	float radius;
	float color; 		// display color
};

// initialize randomly in a box
FUNC_EACH(init,

	p.r = make_xyz_uniform() * W; 	//Random position
	p.radius = rnd_uniform(); 	//Random size (5-7)
	p.density = rnd_normal()+0.5f; 		//Random density (0.5-1.5)

	p.v = make_xyz_uniform()*5;			//Random velocity
)

// linear repulsion + attraction at distance
FUNC_PAIR(pair,
	float ratio = dr / range;
	xyz f = u * (REPULSION_FORCE * (1-ratio) - ATTRACTION_FORCE * ratio);
	addVector(p1.f, f);
	addVector(p2.f, -f);
)


// buoyancy 
FUNC_EACH(buoyancy,
	float volume = p.radius * p.radius * PI;
	p.f.y += (p.density - FLUID_DENSITY) * G * volume;
)

// Euler integration
FUNC_EACH(integrate,
	p.color = xyz_len(p.f) / 50.0f;

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
		fx->runEach(integrate(), A);
		fx->runEach(boundary(), A);

		printf("%.1f ms\n", fx->getTimer());

		// only output to file every 10th step
		if (i % 10 == 0) fx->outputFrame("sample");
	}

	delete fx;
}

