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
        else if (str_eq(key,"energyParticleCount"))         g.energyParticleCount = stoi(val);
        else if (str_eq(key,"nSteps"))                      g.nSteps = stoi(val);
        else if (str_eq(key,"nGenomeInputs"))               g.nGenomeInputs = stoi(val);
        else if (str_eq(key,"interactionRange"))            g.interactionRange = stof(val);
        else if (str_eq(key,"moveFactor"))                  g.moveFactor = stoi(val);
        else if (str_eq(key,"moveCost"))                    g.moveCost = stof(val);
        else if (str_eq(key,"nerveCost"))                   g.nerveCost = stof(val);
        else if (str_eq(key,"genomeCost"))                  g.genomeCost = stof(val);
        else if (str_eq(key,"repulsiveForce"))              g.repulsiveForce = stoi(val);
        else if (str_eq(key,"springForce"))                 g.springForce = stof(val);
        else if (str_eq(key,"groundRepulsiveForce"))        g.groundRepulsiveForce = stof(val);
        else if (str_eq(key,"initialCellEnergy"))           g.initialCellEnergy = stof(val);
        else if (str_eq(key,"minCellEnergy"))               g.minCellEnergy = stof(val);
        else if (str_eq(key,"minDetritusEnergy"))           g.minDetritusEnergy = stof(val);
        else if (str_eq(key,"energyParticleEnergy"))        g.energyParticleEnergy = stof(val);
        else if (str_eq(key,"energyParticleRadius"))        g.energyParticleRadius = stof(val);
        else if (str_eq(key,"cellExistenceThreshold"))      g.cellExistenceThreshold = stof(val);
        else if (str_eq(key,"cellMetabolism"))              g.cellMetabolism = stof(val);
        else if (str_eq(key,"cellDecayRate"))               g.cellDecayRate = stof(val);
        else if (str_eq(key,"fluidDensity"))                g.fluidDensity = stof(val);
        else if (str_eq(key,"gravity"))                     g.gravity = stof(val);
        else if (str_eq(key,"saveFreq"))                    g.saveFreq = stoi(val);
        else if (str_eq(key,"saveIntervalLength"))          g.saveIntervalLength = stoi(val);
        else if (str_eq(key, "saveIntervalDistance"))       g.saveIntervalDistance = stoi(val);
        else if (str_eq(key, "orgInitHealth"))              g.orgInitHealth = stof(val);
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
        stringstream buffer;
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
        /*
        xyz origin = make_xyz_uniform() * int3_to_xyz(g.w);
        origin.y /= 2;
        origin.y += g.w.y / 2;
        */

        //Position in center:
        xyz origin = make_xyz(g.w.x / 2, g.w.y / 2, g.w.z / 2);

        spawnOrganism(
            origin, particleBuffer, p, genome, nerveSys, organisms
            );
    }
    else {
        cerr << "Could not load organism from file: " << path << endl;
    }
}

void loadGenomeAndNerves(string path, Genome *g, NerveSystem *n, int *parent, int *step)
{
    ifstream t(path);
    if (t){
        stringstream buffer;
        buffer << t.rdbuf();

        string org = buffer.str();

        smatch match;

        regex genome_regex("\"genome\":\\{(.*)\\}");
        regex nerve_regex("\"nervesystem\":\\{(.*)\\}");
        regex parent_regex("\"parent\":(.*)");
        regex step_regex("\"step\":(.*)");

        regex_search(org, match, genome_regex);
        string sGenome = match[1];
        *g = Genome(sGenome);

        regex_search(org, match, nerve_regex);
        string sNerveSys = match[1];
        *n = NerveSystem(sNerveSys);

        regex_search(org, match, parent_regex);
        *parent = stoi(match[1]);

        regex_search(org, match, step_regex);
        *step = stoi(match[1]);
    }
    else {
        cerr << "Could not load organism from file: " << path << endl;
    }
}

void logOrgDeath(int organismID, int deathStep, string cause) {
    mkdir("orgDeaths");
    char out_name[256];
    sprintf(out_name, "orgDeaths/org%d.json", organismID);

    FILE *out = fopen(out_name, "wb");
    if (!out) {
        cerr << "Error writing to organism output file " << out_name << endl;
    }

    fprintf(out, "{\n");
    fprintf(out, "\"step\": ");
    fprintf(out, to_string(deathStep).c_str());
    fprintf(out, ",\n\"cause\": \"");
    fprintf(out, cause.c_str());
    fprintf(out, "\"\n}");

    fclose(out);
}

void outputOrganism(Organism *o, int organismID, int creationStep){
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
    fprintf(out, ",\n\"step\": ");
    fprintf(out, to_string(creationStep).c_str());
    fprintf(out, "\n}");

    fclose(out);
}

void dumpCompleteState(Particle *p, int nParticles, int step){
    mkdir("stateDump");
    FILE *out_particles = fopen("stateDump/particles.csv", "wb");
    if (!out_particles) {
        perror("Cannot open file: ");
        QUIT("error opening output file %s\n", "particles.csv");
    }

    FILE *out_step = fopen("stateDump/step.txt", "wb");
    fprintf(out_step, "%i", step);
    fclose(out_step);

    fprintf(out_particles, "particleType,r.x,r.y,r.z,v.x,v.y,v.z,f.x,f.y,f.z,color,radius,alpha,density,energy,energyIn,energyOut,maxEnergy,signal,metabolism,organism,toBuffer,link0,link1,link2,link3,link4,link5,type\n");

    for (int i = 0; i < nParticles; i++) {
        fprintf(
            out_particles,
            "%i,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%i,%d,%i,%i,%i,%i,%i,%i,%i\n",
            p[i].particleType, p[i].r.x, p[i].r.y, p[i].r.z, p[i].v.x, p[i].v.y, p[i].v.z, p[i].f.x, p[i].f.y, p[i].f.z,
            p[i].color, p[i].radius, p[i].alpha, p[i].density, p[i].energy, p[i].energyIn, p[i].energyOut, p[i].maxEnergy,
            p[i].signal, p[i].metabolism, p[i].organism, p[i].toBuffer, p[i].links[0], p[i].links[1], p[i].links[2], p[i].links[3],
            p[i].links[4], p[i].links[5], p[i].type
            );
    }

    fclose(out_particles);
}

int indexOf(string s, vector<string> v){
    ptrdiff_t i = find(v.begin(), v.end(), s) - v.begin();
    if (i >= v.size()) i = -1;
    return (int)i;
}

// Recursive string split
vector<string> split(string s, string delim) {
    //Allocate return vector
    vector<string> items;
    items.reserve(s.size());

    //Find first delimiter position
    size_t delimPos = s.find_first_of(delim);

    //If found, add everything before it as first item
    //and call split recursively with the rest
    if (delimPos != string::npos) {
        //Add first item to vector
        string item = s.substr(0, delimPos);
        items.insert(items.end(), item);
        //Add the rest of the items to vector
        vector<string> rest = split(s.substr(delimPos + delim.length()), delim);
        items.insert(items.end(), rest.begin(), rest.end());
    }
    //If no delimiter is found, this is the last item.
    else {
        items.insert(items.end(), s);
    }
    return items;
}

void loadCompleteState(OrganismMap *organisms, Particle *p, ParticleBuffer *particleBuffer, Fluidix<> *fx, int pSet, int *step, int *currGenomeIndex) {
    ifstream stepStream("stateDump/step.txt");
    if (stepStream) {
        stringstream buffer;
        buffer << stepStream.rdbuf();
        string stepString = buffer.str();
        *step = stoi(stepString);
    }
    else
    {
        cerr << "Could not find file stateDump/step.txt" << endl;
    }

    ifstream particlesFileStream("stateDump/particles.csv");
    std::string line;
    getline(particlesFileStream, line);
    vector<string> headings = split(line, ",");

    int i = 0;
    while (getline(particlesFileStream, line))
    {
        if (i >= fx->getParticleCount(pSet)){
            fx->resizeParticleSet(pSet, i);
            printf("Loading state, increasing particle count to %i", i);
        }
        vector<string> items = split(line, ",");

        p[i].particleType = (ParticleType)stoi(items[indexOf("particleType", headings)]);
        p[i].r.x = stof(items[indexOf("r.x", headings)]);
        p[i].r.y = stof(items[indexOf("r.y", headings)]);
        p[i].r.z = stof(items[indexOf("r.z", headings)]);
        p[i].v.x = stof(items[indexOf("v.x", headings)]);
        p[i].v.y = stof(items[indexOf("v.y", headings)]);
        p[i].v.z = stof(items[indexOf("v.z", headings)]);
        p[i].f.x = stof(items[indexOf("f.x", headings)]);
        p[i].f.y = stof(items[indexOf("f.y", headings)]);
        p[i].f.z = stof(items[indexOf("f.z", headings)]);
        p[i].color = stoi(items[indexOf("color", headings)]);
        p[i].radius = stof(items[indexOf("radius", headings)]);
        p[i].alpha = stof(items[indexOf("alpha", headings)]);
        p[i].density = stof(items[indexOf("density", headings)]);
        p[i].energy = stof(items[indexOf("energy", headings)]);
        p[i].energyIn = stof(items[indexOf("energyIn", headings)]);
        p[i].energyOut = stof(items[indexOf("energyOut", headings)]);
        p[i].maxEnergy = stof(items[indexOf("maxEnergy", headings)]);
        p[i].signal = stof(items[indexOf("signal", headings)]);
        p[i].metabolism = stof(items[indexOf("metabolism", headings)]);
        p[i].organism = stoi(items[indexOf("organism", headings)]);
        p[i].toBuffer = stoi(items[indexOf("toBuffer", headings)]);
        p[i].links[0] = stoi(items[indexOf("link0", headings)]);
        p[i].links[1] = stoi(items[indexOf("link1", headings)]);
        p[i].links[2] = stoi(items[indexOf("link2", headings)]);
        p[i].links[3] = stoi(items[indexOf("link3", headings)]);
        p[i].links[4] = stoi(items[indexOf("link4", headings)]);
        p[i].links[5] = stoi(items[indexOf("link5", headings)]);
        p[i].type = (CellType)stoi(items[indexOf("type", headings)]);

        //If cell, build up organism map
        if (p[i].particleType == Cell && p[i].organism != -1) {
            Organism o;
            try {
                o = organisms->at(p[i].organism);
            }
            catch (const std::out_of_range& e) {
                Genome genome;
                NerveSystem nervousSystem;
                int parent, birthStep;

                char path[256];
                sprintf(path, "organisms/org%d.json", p[i].organism);
                loadGenomeAndNerves(path, &genome, &nervousSystem, &parent, &birthStep);

                o.genome = genome;
                o.nerveSystem = nervousSystem;
                o.parent = parent;
                o.health = g.orgInitHealth - g.dt * (*step - birthStep);

                if (*currGenomeIndex < p[i].organism)
                    *currGenomeIndex = p[i].organism;
            }
            o.cells.push_back(i);
            int orgID = p[i].organism;
            organisms->operator[](orgID) = o;
        }
        //If instead buffer, build upp the buffer queue
        else if (p[i].particleType == Buffer){
            particleBuffer->push(i);
        }
        i++;
    }
    fx->resizeParticleSet(pSet, i);
    g.nParticles = i;
}

void outputParticles(Particle *p, int nParticles, int step) {
    char out_name[256];
    sprintf(out_name, "output/frame%d.csv", step);

    mkdir("output");
    FILE *out = fopen(out_name, "wb");
    if (!out) {
        perror("Cannot open file: ");
        QUIT("error opening output file %s\n", out_name);
    }

    fprintf(out, "pt,ct,o,e,r,x,y,z\n");
    for (int i = 0; i < nParticles; i++) {
        if (p[i].particleType == Cell || p[i].particleType == Detritus) {
            fprintf(out,"%i,%i,%i,%f,%f,%f,%f,%f\n",
                (int)p[i].particleType,
                (int)p[i].type,
                p[i].organism,
                p[i].energy,
                p[i].radius,
                p[i].r.x, p[i].r.y, p[i].r.z
            );
        }
    }
    fclose(out);
}

void saveTerrain(Particle *p, vector<tuple<int, int>> links, int nParticles) {
    mkdir("terrain");
    FILE *outMesh = fopen("terrain/mesh.csv", "wb");
    if (!outMesh) {
        perror("Cannot open file: ");
        QUIT("error opening output file: terrain/mesh.csv\n");
    }

    fprintf(outMesh, "x,y,z\n");
    for (int i = 0; i < nParticles; i++) {
        fprintf(outMesh, "%f,%f,%f\n",
            p[i].r.x, p[i].r.y, p[i].r.z
        );
    }
    fclose(outMesh);

    FILE *outMeshLinks = fopen("terrain/links.csv", "wb");
    if (!outMeshLinks) {
        perror("Cannot open file: ");
        QUIT("error opening output file: terrain/links.csv\n");
    }

    fprintf(outMeshLinks, "i,o\n");
    for (int i = 0; i < links.size(); i++) {
        fprintf(outMeshLinks, "%i,%i\n",
            get<0>(links[i]),
            get<1>(links[i])
        );
    }
    fclose(outMeshLinks);
}

#endif