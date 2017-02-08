#ifndef IO_H
#define IO_H

#include <fstream>
#include <sstream>
#include <algorithm>
#include <concurrent_queue.h>
#include "../lib/structures.h"

using namespace std;

//Define spawnOrganism function
int spawnOrganism(
    xyz origin, ParticleBuffer *particleBuffer, Particle *p,
    Genome genome, NerveSystem nerveSys, OrganismMap *organisms
);

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

Global loadConfig(string path) {
    ifstream file;
    file.open(path.c_str());
    if (!file){
        printf("File %s could not be found!\n", path.c_str());
        exit(-1);
    }
    else
       printf("Loaded file: %s:\n", path.c_str());

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

void loadOrg(
    string path,
    ParticleBuffer *particleBuffer,
    Particle *p,
    OrganismMap *organisms)
{
    ifstream t(path);
    if (t){

        std::stringstream buffer;
        buffer << t.rdbuf();

        string org = buffer.str();

        smatch match;

        regex genome_regex("\"genome\":\\{(.*)\\}");
        regex nerve_regex("\"nervesystem\":\\{(.*)\\}");
        //regex parent_regex("\"parent\":(.*)");

        regex_search(org, match, genome_regex);
        string sGenome = match[1];
        Genome genome(sGenome);

        regex_search(org, match, nerve_regex);
        string sNerveSys = match[1];
        NerveSystem nerveSys(sNerveSys);

        //regex_search(org, match, parent_regex);
        //int parent = stoi(match[1]);

        //Set random position:
        xyz origin = make_xyz_uniform() * int3_to_xyz(g.w);
        origin.y /= 2;
        origin.y += g.w.y / 2;

        spawnOrganism(
            origin, particleBuffer, p, genome, nerveSys, organisms
            );
    }
    else {
        cerr << "Could not load organism from file: " << path << endl;
    }
}

void outputOrganism(Organism *o, int organismID){
    mkdir("organisms");
    char out_name[256];
    sprintf(out_name, "organisms/org%d.json", organismID);

    FILE *out = fopen(out_name, "wb");
    if (!out) {
        cerr << "Error writing to organism output file " << out_name << endl;
    }

    fprintf(out, "{\n");
    fprintf(out, o->genome.toJSON().c_str());
    fprintf(out, ",\n");
    fprintf(out, o->nerveSystem.toJSON().c_str());
    fprintf(out, ",\n\"parent\": ");
    fprintf(out, to_string(o->parent).c_str());
    fprintf(out, "\n}");

    fclose(out);
}


void outputParticles(Particle *p, int nParticles, int step) {
    char out_name[256];
    sprintf(out_name, "output/frame%d.json", step);

    mkdir("output");
    FILE *out = fopen(out_name, "wb");
    if (!out) {
        perror("Cannot open file: ");
        QUIT("error opening output file %s\n", out_name);
    }

    fprintf(out, "{\"Items\":[\n");
    bool first = true;
    for (int i = 0; i < nParticles; i++) {
        if (p[i].particleType == Cell || p[i].particleType == Pellet) {
            if (first){
                first = false;
            }
            else {
                fprintf(out, ",");
            }
            fprintf(out,
                //"pt:%i, ct:%i, o:%i, x:%f, y:%f, z:%f",
                "{\"pt\":%i,\"ct\":%i,\"o\":%i,\"e\":%f,\"x\":%f,\"y\":%f,\"z\":%f}\n",
                (int)p[i].particleType,
                (int)p[i].type,
                p[i].organism,
                p[i].energy,
                p[i].r.x, p[i].r.y, p[i].r.z
                );
        }
    }
    fprintf(out, "]}");
    fclose(out);
}

#endif