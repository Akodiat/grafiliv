#include "fluidix.h"
#include "../lib/genome.h"

#define L 8

struct Global {
	float grid_pack;
	int3 grid_num;
	xyz grid_origin;
} g;

struct Particle {
	xyz r, v, f;
	float color;
	float alpha;
};

FUNC_EACH(put_in_grid,
	int3 ijk;
	ijk.x = p_index % g.grid_num.x;
	ijk.y = (p_index / g.grid_num.x) % g.grid_num.y;
	ijk.z = p_index / (g.grid_num.x * g.grid_num.y);

	p.r = g.grid_origin + g.grid_pack * int3_to_xyz(ijk);
)

void createParticles(Fluidix<> *fx, xyz origin, int side){
	xyz box_min = origin;
	xyz box_max = origin + make_xyz(side, side, side);

	g.grid_pack = 1.0f; // distance between adjacent particles

	xyz box_size = box_max - box_min;

	g.grid_num.x = (int)roundf(box_size.x / g.grid_pack);
	g.grid_num.y = (int)roundf(box_size.y / g.grid_pack);
	g.grid_num.z = (int)roundf(box_size.z / g.grid_pack);

	g.grid_origin = (box_max + box_min - g.grid_pack * (int3_to_xyz(g.grid_num) - 1)) / 2;

	int setA = fx->createParticleSet(g.grid_num.x * g.grid_num.y * g.grid_num.z);
	fx->runEach(put_in_grid(), setA);

	int inputs = 4;
	int outputs = 2;
	Genome g(inputs, outputs);
	g.mutate();
	g.printMathematica();

	Particle *p = fx->getParticleArray(setA);
	int n = fx->getParticleCount(setA);

	for(int step=0; step<20; step++) {
		g.mutate();
		for (int i = 0; i < n; i++) {
			xyz dr = p[i].r - make_xyz(L/2, L/2, L/2);
	
			vector<float> input;
			input.push_back(dr.x);
			input.push_back(dr.y);
			input.push_back(dr.z);
			input.push_back(xyz_len(dr));
	
			vector<float> output = g.getOutput(input);
	
			printf("input: ");  for(float i : input)  printf("%.2f ",i); printf("\t");
			printf("output: "); for(float o : output) printf("%.2f ",o); printf("\n");
	
			p[i].color = output[0];
			p[i].alpha = output[1]>0 ? 1.0f : 0.0f;
	
			// apply changes before any operation
			fx->applyParticleArray(setA);
		}
		g.printMathematica();
		g.printGenome();
		fx->outputFrame("output");
	}
}

int main() {
	Fluidix<> *fx = new Fluidix<>(&g);
	
	createParticles(fx, make_xyz(0, 0, 0), L);


	delete fx;
}
