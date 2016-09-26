#include "fluidix.h"
#include "../lib/genome.h"
#include <map>

#define W 100 // simulation size
#define SUNLIGHT_STRENGTH 0.1f
#define MAX_CELL_ENERGY 5.0f
#define DT 0.01f
#define RANGE 3.0f

#define REPULSION_FORCE 150
#define ATTRACTION_FORCE 40

#define FLUID_DENSITY 1.0f
#define G -9.81f

using namespace std;

__host__ __device__ bool operator==(xyz a, xyz b) {
	return (
		a.x == b.x && 
		a.y == b.y &&
		a.z == b.z
	);
}

struct Global {
} g;

struct Particle {
	xyz r, v, f;
	float energy;
	float color;
	float alpha;
	float radius;
	float density;
	bool remove;
	bool divide;
	int linkSet;
	int origin;
};


FUNC_EACH(init,
	p.r = make_xyz_uniform() * W;
	p.energy = 1000.0f;
	p.alpha  = 1.0f;
	p.density = FLUID_DENSITY;
	p.radius = 1.0f;
	p.remove = false;
	p.divide = false;
	p.linkSet = -1;
	p.origin = p_index;
)

// pair interaction
FUNC_PAIR(particlePair,
	float ratio = dr / range;
	float attraction = 0.0f;
	if(p1.linkSet == p2.linkSet)
		attraction = ATTRACTION_FORCE * ratio;
	xyz f = u * (REPULSION_FORCE * (1-ratio) - attraction);
	addVector(p1.f, f);
	addVector(p2.f, -f);
)

// buoyancy 
FUNC_EACH(buoyancy,
	float volume = p.radius * p.radius * PI;
	p.f.y += (p.density - FLUID_DENSITY) * G * volume;
)

FUNC_EACH(integrate,
	p.v += p.f * DT;
	p.r += p.v * DT;
	p.f = make_xyz(0, 0, 0);
)

FUNC_EACH(photosynthesis,
	p.energy += SUNLIGHT_STRENGTH;
)

FUNC_EACH(friction,
	p.v *= 0.9f;
)

FUNC_EACH(checkDivide,
	if (p.energy > MAX_CELL_ENERGY && rnd_uniform() < 0.1f)
		p.divide = true;
	else if(p.energy <= 0)
		p.remove = true;	
)

FUNC_EACH(boundary,
	if (p.r.x < 0) { p.v.x = 0.9f * (0 - p.r.x) / DT; p.r.x = 0; }
	if (p.r.x > W) { p.v.x = 0.9f * (W - p.r.x) / DT; p.r.x = W; }
	if (p.r.y < 0) { p.v.y = 0.9f * (0 - p.r.y) / DT; p.r.y = 0; }
	if (p.r.y > W) { p.v.y = 0.9f * (W - p.r.y) / DT; p.r.y = W; }
	if (p.r.z < 0) { p.v.z = 0.9f * (0 - p.r.z) / DT; p.r.z = 0; }
	if (p.r.z > W) { p.v.z = 0.9f * (W - p.r.z) / DT; p.r.z = W; }
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
		fx->runEach(friction(), A);
		fx->runEach(checkDivide(), A);
		fx->runEach(boundary(), A);

		if (step % 1 == 0) {
			//printf("step %d\n", step);
			fx->outputFrame("output");
		}
	}

	delete fx;
}

