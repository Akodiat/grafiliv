#ifndef LOAD_CONFIG_H
#define LOAD_CONFIG_H

#include <fstream>
#include <algorithm>

using namespace std;

struct Global {
    int nEggs;

    float dt; // integration time-step
    float pelletLifetime;
    int3 w;
    int nParticles;
    int nInitialOrganisms;
    int bufferSize;
    int nSteps;
    int3 initialOrganismDimensions;
    int nGenomeInputs; //Inputs x,y,z,d
    float interactionRange;
    int moveFactor;
    int repulsiveForce;
    float springForce; // spring constant
    float groundRepulsiveForce; // repulsive wall force
    float initialCellEnergy;
    float minCellEnergy;
    float minPelletEnergy;
    float energyParticleEnergy;
    float energyParticleRadius;
    float cellExistenceThreshold;
    float cellMetabolism;
    float cellDecayRate;
    float fluidDensity;
    float gravity;
} g;

// Convert string to int3
// Requires form "x,y,z"
int3 stoi3(string s) {
    size_t pos1 = s.find(",");
    size_t pos2 = s.rfind(",");

    int x = stoi(s.substr(0, pos1));
    int y = stoi(s.substr(pos1 + 1, pos2 - 1));
    int z = stoi(s.substr(pos2 + 1, s.length()));

    return make_int3(x, y, z);
}

bool str_eq(string a, string b){
    return !a.compare(b);
}

Global load(string path) {
    ifstream file;
    file.open(path.c_str());
    if (!file){
        printf("Config file %s could not be found!\n", path.c_str());
        exit(-1);
    }
    else
       printf("Loaded config file: %s:\n", path.c_str());

    Global g = Global();

    string line;
    while (getline(file, line))
    {
        //Disregard everything after the semicolon:
        size_t commentPos = line.find(";");
        string s = line.substr(0, commentPos);

        //Remove spaces:
        s.erase(remove_if(s.begin(), s.end(), ::isspace), s.end());

        //Disregard empty lines:
        if (s.length() == 0) continue;

        size_t pos = s.find("=");

        string key = s.substr(0, pos);
        string val = s.substr(pos + 1, s.length());

        printf("parameter: %s, value: %s\n", key.c_str(), val.c_str());

        if      (str_eq(key,"dt"))                          g.dt = stof(val);
        else if (str_eq(key,"pelletLifetime"))              g.pelletLifetime = stof(val);  
        else if (str_eq(key,"w"))                           g.w = stoi3(val);  
        else if (str_eq(key,"nParticles"))                  g.nParticles = stoi(val);  
        else if (str_eq(key,"nInitialOrganisms"))           g.nInitialOrganisms = stoi(val);  
        else if (str_eq(key,"bufferSize"))                  g.bufferSize = stoi(val);  
        else if (str_eq(key,"nSteps"))                      g.nSteps = stoi(val);  
        else if (str_eq(key,"initialOrganismDimensions"))   g.initialOrganismDimensions = stoi3(val);  
        else if (str_eq(key,"nGenomeInputs"))               g.nGenomeInputs = stoi(val);  
        else if (str_eq(key,"interactionRange"))            g.interactionRange = stof(val);  
        else if (str_eq(key,"moveFactor"))                  g.moveFactor = stoi(val);  
        else if (str_eq(key,"repulsiveForce"))              g.repulsiveForce = stof(val);  
        else if (str_eq(key,"springForce"))                 g.springForce = stof(val);  
        else if (str_eq(key,"groundRepulsiveForce"))        g.groundRepulsiveForce = stof(val);  
        else if (str_eq(key,"initialCellEnergy"))           g.initialCellEnergy = stof(val);  
        else if (str_eq(key,"minCellEnergy"))               g.minCellEnergy = stof(val);  
        else if (str_eq(key,"minPelletEnergy"))             g.minPelletEnergy = stof(val);  
        else if (str_eq(key,"energyParticleEnergy"))        g.energyParticleEnergy = stof(val);
        else if (str_eq(key,"energyParticleRadius"))        g.energyParticleRadius = stof(val);
        else if (str_eq(key,"cellExistenceThreshold"))      g.cellExistenceThreshold = stof(val);  
        else if (str_eq(key,"cellMetabolism"))              g.cellMetabolism = stof(val);  
        else if (str_eq(key,"cellDecayRate"))               g.cellDecayRate = stof(val);  
        else if (str_eq(key,"fluidDensity"))                g.fluidDensity = stof(val);  
        else if (str_eq(key,"gravity"))                     g.gravity = stof(val);  
    }

    file.close();

    return g;
}

#endif