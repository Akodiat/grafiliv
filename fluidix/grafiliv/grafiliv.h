#ifndef GRAFILIV_H
#define GRAFILIV_H

// Initialize particle as cell
void setDefaultCellValues(Particle *cell);

// Given a phenotype network output, apply it to the cell
bool applyPhenotype(vector<float> output, Particle *cell);

// Remove cell links from and to cell
void disconnectCell(Particle *p, int cell, int code);

void emptyCellPos(Particle *p, int cell);
void deadCellPos(Particle *p, int cell) ;

// Helper function to get the 1-dimensional index
// given x,y,z and box size br
int getIdxFromCoord(int x, int y, int z, int3 br);

// Create cells of an organism given a genome and a nervous system
pair<int, vector<int>> createCellsFromGenotype(
    xyz origin, ParticleBuffer *particleBuffer,
    Particle *p, Genome *genome, NerveSystem *nerveSys, 
	OrganismMap *organisms
);

//Initialize new organism (without parent)
int spawnOrganism(
    xyz origin, ParticleBuffer *particleBuffer,
    Particle *p, Genome genome, NerveSystem nerveSys, OrganismMap *organisms);
	
	// Initialize new organism from parent
int spawnOrganism(
    xyz origin, ParticleBuffer *particleBuffer,
    Particle *p, int parent, OrganismMap *organisms);
	
// Initialize a random organism
int initializeOrganism(ParticleBuffer *particleBuffer, Particle *p, OrganismMap *organisms);

Matrix3 getTransform(xyz front, xyz right, xyz up, xyz back, xyz left, xyz down);

#endif