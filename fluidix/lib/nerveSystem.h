#ifndef NERVE_SYSTEM_H
#define NERVE_SYSTEM_H

#define MUTATION_PROB 0.01f

#include <random>
#include <time.h>
#include <regex>
#include <map>

using namespace std;

//default_random_engine rndGen(time(0));

class NerveSystem {

public:
    enum NodeType { Input, Hidden, Output };

    struct Node {
        NodeType type;
        float preVal;
        float postVal;

        Node(NodeType t) :
            type(t), preVal(0.0f), postVal(0.0f)
        {}

    Node(){}

        float activationFunction(float x) {
        return x / (1.0f + abs(x));
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
    NerveSystem() {}

    // Create nervous system, specifying number of input and output nodes
    NerveSystem(int nOut) {
        nOutputs = nOut;
    nInputs = 0;
    nextNodeId = 0;

        // Add output nodes
        for(int i=0; i<nOut; i++) {
            nodes.emplace(nextNodeId++, Node(Output));
        }
    }

    NerveSystem(string json) {
        // Extraction of several sub-matches
        regex vertices_regex("\"vertices\":\\[(.*?)\\]");
        regex links_regex("\"links\":\\[(.*?)\\]");
        regex vertex_regex("\"i\":(\\d+),\"type\":(\\d+)");
        regex link_regex("\"i\":(\\d+),\"o\":(\\d+),\"w\":(-?\\d+(.\\d+)?)");

        smatch match;

        regex_search(json, match, vertices_regex);
        string vertices = match[1];

        regex_search(json, match, links_regex);
        string links = match[1];

        while (regex_search(vertices, match, vertex_regex)){
            ssub_match idMatch = match[1];
            ssub_match typeMatch = match[2];

            NodeType type = (NodeType)stoi(typeMatch.str());
            int id = stoi(idMatch.str());

            Node node(type);
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

        nInputs = nOutputs = 0;
        for(int i=0; i<nodes.size(); i++) {
            switch (nodes[i].type) {
            case Input:     nInputs++; break;
            case Output:    nOutputs++; break;
            default:        break;
            }
        }
    }

    void updateInputs(int nIn) {
        normal_distribution<float> rndNormal(0.0f, 1.0f);
        while (nInputs < nIn) {
            int in = nextNodeId;
            nodes.emplace(nextNodeId++, Node(Input));
            nInputs++;
            for(int i=0; i<nodes.size(); i++) {
                if (nodes[i].type == Output) {
                    float weight = rndNormal(rndGen);
                    connections.push_back(Connection(in, i, weight, true));
                }
            }
        }

    }

    //Get size of genome
    int getSize(){
        return connections.size() +
            nodes.size() - nInputs - nOutputs;
    }

    // Input data into the network and get output
    vector<float> getOutput(vector<float> input) {
        // Clear previous values
//      for (int i = 0; i<nodes.size(); i++) {
//          nodes[i].preVal = nodes[i].postVal = 0.0f;
//      }

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
                "{\"i\":"               +
                to_string(i)           +
                ",\"type\":"            +
                to_string(nodes[i].type)    +
                "}";
        }

        first = true;
        for(int i=0; i<connections.size(); i++) {
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
            string("\"nervesystem\":{") +
            string("\"vertices\":[") + vertices + string("],") +
            string("\"links\":[") + links + string("]}");
    }

    private:
        vector<Connection>          connections;
        map<int, Node>  nodes;
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
            for (int i = 0; i<connections.size(); i++)
            if (rndUniform(rndGen) < MUTATION_PROB)
                connections[i].weight *= rndNormal(rndGen);
        }

        // Mutate structure by adding connection
        void mutateAddConnection() {
            // If there are connection left to add
            if (connections.size() < nodes.size()*nodes.size()) {
                vector<int> sources, recievers;
                for (auto it : nodes) {
                    int id = it.first;
                    Node n = it.second;
                    {
                        switch (n.type) {
                        case Input:
                            sources.push_back(id);
                            break;
                        case Output:
                            recievers.push_back(id);
                            break;
                        case Hidden:
                            sources.push_back(id);
                            recievers.push_back(id);
                            break;
                        }
                    }
                }

                if (sources.size() == 0 || recievers.size() == 0)
                    return;

                uniform_int_distribution<int> rndSource(0, sources.size() - 1);
                uniform_int_distribution<int> rndReciever(0, recievers.size() - 1);

                int in = sources[rndSource(rndGen)];
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


        // Propagate values through network one step
        void computeOneCycle(){
            for (int i = 0; i<connections.size(); i++) {
                if (connections[i].expressed) {
                    nodes[connections[i].out].preVal +=
                        nodes[connections[i].in].postVal *
                        connections[i].weight;
                }
            }
            for (auto it : nodes) {
                Node node = it.second;
                if (node.type != Input) {
                    node.postVal = node.activationFunction(node.preVal);
                }
                node.preVal = 0.0f;
            }
        }
};

int NerveSystem::currInnovNumber = 0;

#endif
