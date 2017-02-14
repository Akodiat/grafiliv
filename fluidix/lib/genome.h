#ifndef GENOME_H
#define GENOME_H

#define MUTATION_PROB 0.01f

#include <random>
#include <time.h>
#include <regex>
#include <map>

using namespace std;

default_random_engine rndGen(time(0));

class Genome {
public:
    enum NodeType { Input, Hidden, Output };
    enum ActivationFunction {
        Sine, Abs, Id, Mod, Gaus,
        N_ACTIVATION_FUNCTIONS
    };
    struct Node {
        NodeType type;
        ActivationFunction f;
        float preVal;
        float postVal;

        Node(NodeType t) :
            type(t), preVal(0.0f), postVal(0.0f)
        {
            uniform_int_distribution<int> rndInt(0, N_ACTIVATION_FUNCTIONS - 1);
            f = ActivationFunction(rndInt(rndGen));
        }

        Node(NodeType t, ActivationFunction func)
            //: type(t), preVal(0.0f), postVal(0.0f), f(func)
        {
            type = t;
            f = func;
            preVal = 0.0f;
            postVal = 0.0f;
        }

        Node(){}

        float activationFunction(float x) {
            switch (f) {
            case Sine:  return sin(x);
            case Abs:   return clamp(1.0f - abs(x), -1, 1);
            case Gaus:  return exp(-(x*x) / 2);
            case Id:    return clamp(x, -1, 1);
            case Mod:   return fmod(x, 1);
            default:
                printf("Activation function (%d) not set", f);
                return NULL;
            }
        }
    };

    struct Connection {
        int in;
        int out;
        float weight;
        bool expressed;
        int innovNumber;
        Connection(int i, int o, float w, bool e) :
            in(i), out(o), weight(w), expressed(e),
            innovNumber(currInnovNumber++) {}
    };
    // Default constructor, needed to initialize array
    Genome() {}

    // Create genome, specifying number of input and output nodes
    Genome(int nIn, int nOut, int3 boundRadius) {
        normal_distribution<float> rndNormal(0.0f, 1.0f);

        boundingRadius = boundRadius;
        nInputs = nIn;
        nOutputs = nOut;
        // Add input nodes
        for (int i = 0; i<nIn; i++) {
            nodes.emplace(nextNodeId++, Node(Input));
        }
        // Add output nodes and connect them to each input node
        for (int i = 0; i<nOut; i++) {
            nodes.emplace(nextNodeId++, Node(Output));
            int out = nodes.size() - 1;
            for (int j = 0; j<nIn; j++) {
                int in = j;
                float weight = rndNormal(rndGen);
                connections.push_back(Connection(in, out, weight, true));
            }
        }
    }

    Genome(string json) {
        // Extraction of several sub-matches
        regex vertices_regex("\"vertices\":\\[(.*?)\\]");
        regex links_regex("\"links\":\\[(.*?)\\]");
        regex radius_regex("\"radius\":\\[(.*?),(.*?),(.*?)\\]");
        regex vertex_regex("\"i\":(\\d+),\"type\":(\\d+),\"f\":(\\d+)");
        regex link_regex("\"i\":(\\d+),\"o\":(\\d+),\"w\":(-?\\d+(.\\d+)?)");

        smatch match;

        regex_search(json, match, vertices_regex);
        string vertices = match[1];

        regex_search(json, match, links_regex);
        string links = match[1];

        while (regex_search(vertices, match, vertex_regex)){
            ssub_match idMatch = match[1];
            ssub_match typeMatch = match[2];
            ssub_match fMatch = match[3];

            NodeType type = (NodeType)stoi(typeMatch.str());
            int id = stoi(idMatch.str());
            ActivationFunction f = (ActivationFunction)stoi(fMatch.str());

            Node node(type, f);
            nodes.emplace(id, node);

            vertices = match.suffix().str();
        }

        while (regex_search(links, match, link_regex)){
            ssub_match in       = match[1];
            ssub_match out      = match[2];
            ssub_match weight   = match[3];

            connections.push_back(Connection(
                stoi(in.str()),
                stoi(out.str()),
                stof(weight.str()),
                true));

            links = match.suffix().str();
        }

        if (regex_search(json, match, radius_regex)){
            ssub_match x = match[1];
            ssub_match y = match[2];
            ssub_match z = match[3];
            boundingRadius = make_int3(
                stoi(x.str()),
                stoi(y.str()),
                stoi(z.str())
                );
        }

        nInputs = nOutputs = 0;
        for(int i=0; i<nodes.size(); i++) {
            switch (nodes[i].type) {
            case Input:     nInputs++; break;
            case Output:    nOutputs++; break;
            default:        break;
            }
        }
    }

    int3 getBoundingRadius() {
        return boundingRadius;
    }
    int getMaxCellsReq() {
        int x = 2 * boundingRadius.x + 1;
        int y = 2 * boundingRadius.y + 1;
        int z = 2 * boundingRadius.z + 1;
        return x*y*z;
    }

    //Get size of genome
    int getSize(){
        return connections.size() +
            nodes.size() - nInputs - nOutputs;
    }

    // Input data into the network and get output
    vector<float> getOutput(vector<float> input) {
        // Clear previous values
        for (int i = 0; i<nodes.size(); i++) {
            nodes[i].preVal = nodes[i].postVal = 0.0f;
        }

        int iInput = 0;
        for (int i = 0; i<nodes.size(); i++) {
            if (nodes[i].type == Input){
                // It is possible for a sensor cell to die, in that case
                // its signal value will be 0.0f
                try {
                    nodes[i].postVal = input.at(iInput);
                }
                catch (const std::out_of_range& oor){
                    nodes[i].postVal = 0.0f;
                }
                iInput++;
            }
        }
        for (int i = 0; i<nActivationCycles; i++) {
            computeOneCycle();
        }
        vector<float> output;
        for(int i=0; i<nodes.size(); i++) {
            if (nodes[i].type == Output){
                //TODO: should we use activationFunction here?
                output.push_back(nodes[i].activationFunction(nodes[i].postVal));
            }
        }
        return output;
    }

    // Mutate network
    void mutate(){
        normal_distribution<float> rndNormal(0.0f, 1.0f);
        mutateConnections();
        mutateBoundingRadius();
        if (rndNormal(rndGen) < MUTATION_PROB)
            mutateAddConnection();
        if (rndNormal(rndGen) < MUTATION_PROB)
            mutateAddNode();
        if (rndNormal(rndGen) < MUTATION_PROB)
            mutateRemoveConnection();
    }

    string toJSON() {
        string vertices = "";
        string links = "";

        bool first = true;
        for(int i = 0; i<nodes.size(); i++) {
            if (first) first = false;
            else vertices += ",";

            vertices +=
                "{\"i\":" +
                to_string(i) +
                ",\"type\":" +
                to_string(nodes[i].type) +
                ",\"f\":" +
                to_string(nodes[i].f) +
                "}";
        }

        first = true;
        for (int i = 0; i<connections.size(); i++) {
            if (connections[i].expressed) {
                if (first) first = false;
                else {
                    links += ",";
                }

                links +=
                    "{\"i\":" + to_string(connections[i].in) +
                    ",\"o\":" + to_string(connections[i].out) +
                    ",\"w\":" + to_string(connections[i].weight) +
                    "}";
            }
        }

        return
            string("\"genome\":{") +
            string("\"vertices\":[") + vertices + string("],") +
            string("\"links\":[") + links + string("],") +
            string("\"radius\":[") +
            (
            to_string(boundingRadius.x) + "," +
            to_string(boundingRadius.x) + "," +
            to_string(boundingRadius.x) + "]}"
            );
    }

    void printMathematica() {
        string vertices = "";
        string links = "";
        string vertexStyle = "";
        string edgeWeights = "";
        for (int i = 0; i<nodes.size(); i++) {
            vertices += to_string(i) +
                (i == nodes.size() - 1 ? "" : ", ");
            switch (nodes[i].type) {
            case Input: vertexStyle += (to_string(i) + "->Blue, "); break;
            case Output: vertexStyle += (to_string(i) + "->Red, ");  break;
            }
        }
        for (int i = 0; i<connections.size(); i++) {
            links +=
                to_string(connections[i].in) + "->" +
                to_string(connections[i].out) +
                (i == connections.size() - 1 ? "" : ", ");
            edgeWeights += to_string(connections[i].weight) +
                (i == connections.size() - 1 ? "" : ", ");
        }
        cout << "Graph[{" << vertices << "},{" << links << "}, " <<
            "VertexStyle->{" << vertexStyle << "}, " <<
            "EdgeWeight->{" << edgeWeights << "}, " <<
            "VertexLabels->\"Name\"]" << endl;
    }

private:
    vector<Connection>          connections;
    map<int,Node>               nodes;
    int3                        boundingRadius;
    static int currInnovNumber;
    int nextNodeId;

    int nInputs;
    int nOutputs;

    // Number of times to propagate values
    int nActivationCycles = 10;

    // Mutate each connection with small pertubation
    void mutateConnections() {
        uniform_real_distribution<float> rndUniform(0.0f, 1.0f);
        normal_distribution<float> rndNormal(0.0f, 1.0f);
        for(int i=0; i<connections.size(); i++)
            if (rndUniform(rndGen) < MUTATION_PROB)
                connections[i].weight *= rndNormal(rndGen);
    }

    // Mutate structure by adding connection
    void mutateAddConnection() {
        // If there are connection left to add
        if (connections.size() < nodes.size()*nodes.size()) {
            vector<int> sources, recievers;
            for(int i=0; i<nodes.size(); i++) {
                {
                    switch (nodes[i].type) {
                    case Input:
                        sources.push_back(i);
                        break;
                    case Output:
                        recievers.push_back(i);
                        break;
                    case Hidden:
                        sources.push_back(i);
                        recievers.push_back(i);
                        break;
                    }
                }
            }
            if (sources.size() == 0 || recievers.size() == 0)
                    return;

            uniform_int_distribution<int> rndSource(0, sources.size() - 1);
            uniform_int_distribution<int> rndReciever(0, recievers.size() - 1);

            int in  = sources[rndSource(rndGen)];
            int out = recievers[rndReciever(rndGen)];

            int timeoutTries = sources.size() * recievers.size();

            // Make sure the connection doesn't already exist
            for (int i = 0; i < connections.size(); i++) {
                if (connections[i].in == in && connections[i].out == out) {
                    // For n nodes, give up after n^2 tries
                    if (!timeoutTries--) return;

                    int in = rndSource(rndGen);
                    int out = rndReciever(rndGen);

                    i = 0;
                }
            }
            normal_distribution<float> rndNormal(0.0f, 1.0f);
            float weight = rndNormal(rndGen);

            // Add new connection (with random weight) to list
            connections.push_back(Connection(in, out, weight, true));
        }
    }

    // Mutate structure by adding node
    void mutateAddNode() {
        if (connections.size() > 0) {
            // Select and disable random connection
            uniform_int_distribution<int> rndInt(0, connections.size() - 1);
            int iRndCon = rndInt(rndGen);
            connections[iRndCon].expressed = false;

            nodes.emplace(nextNodeId++, Node(Hidden)); // Add new node to list
            int id = nextNodeId; // Get new node's index

            // Connect new node between previously connected nodes
            connections.push_back(Connection(
                connections[iRndCon].in,
                id,
                1.0f,
                true
                ));
            connections.push_back(Connection(
                id,
                connections[iRndCon].out,
                connections[iRndCon].weight,
                true
                ));
        }
    }

    // Mutate structure by removing connection
    void mutateRemoveConnection() {
        // Select and remove a random node
        if (connections.size() > 0) {
            uniform_int_distribution<int> rndInt(0, connections.size() - 1);
            int connection = rndInt(rndGen);
            connections.erase(connections.begin() + connection);

            for (auto it = nodes.cbegin(); it != nodes.cend();){
                int id = (*it).first;
                Node node = (*it).second;
                if (node.type == Hidden) {
                    //Check for remaining connections to the hidden node;
                    for (int j = 0; j < connections.size(); j++) {
                        if (connections[j].in == id || connections[j].out == id){
                            ++it;
                            continue;
                        }
                    }
                    // If there are no remaining connections to the hidden node
                    // it can be removed:
                    it = nodes.erase(it);
                }
            }
        }
    }

    void mutateBoundingRadius() {
        uniform_real_distribution<float> rndUniform(0.0f, 1.0f);
        //Divide by six to make prob of mutation in any direction
        //into MUTATION_PROB
        if (rndUniform(rndGen) < MUTATION_PROB / 6)
            boundingRadius.x++;
        if (rndUniform(rndGen) < MUTATION_PROB / 6)
            boundingRadius.x = max(boundingRadius.x - 1, 0);
        if (rndUniform(rndGen) < MUTATION_PROB / 6)
            boundingRadius.y++;
        if (rndUniform(rndGen) < MUTATION_PROB / 6)
            boundingRadius.y = max(boundingRadius.y - 1, 0);
        if (rndUniform(rndGen) < MUTATION_PROB / 6)
            boundingRadius.z++;
        if (rndUniform(rndGen) < MUTATION_PROB / 6)
            boundingRadius.z = max(boundingRadius.z - 1, 0);
    }

    // Propagate values through network one step
    void computeOneCycle(){
        for (int i = 0; i<connections.size(); i++) {
            if (connections[i].expressed) {
                nodes[connections[i].out].preVal +=
                    nodes[connections[i].in].postVal *
                    connections[i].weight;
            }
        }
        for(int i=0; i<nodes.size(); i++) {
            if (nodes[i].type != Input) {
                nodes[i].postVal = nodes[i].activationFunction(nodes[i].preVal);
            }
            nodes[i].preVal = 0.0f;
        }
    }
};

int Genome::currInnovNumber = 0;

#endif
