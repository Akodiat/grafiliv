#ifndef NERVE_SYSTEM_H
#define NERVE_SYSTEM_H

#define MUTATION_PROB 0.01f

#include <random>
#include <time.h>

using namespace std;

//default_random_engine rndGen(time(0));

class NerveSystem {
	enum NodeType {Input, Hidden, Output};
	struct Node {
		NodeType type;
		float preVal;
		float postVal;

		Node(NodeType t) :
			type(t), preVal(0.0f), postVal(0.0f)
		{}

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
			in(i), out(o), weight(w), 	expressed(e),
			innovNumber(currInnovNumber++) {}
	};

private:
	vector<Connection>	connections;
	vector<Node> 		nodes;
	static int currInnovNumber;

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
		uniform_real_distribution<float> rndUniform(0.0f, 1.0f);
		normal_distribution<float> rndNormal(0.0f, 1.0f);

		// If there are connection left to add
		if (connections.size() < nodes.size()*nodes.size()) {
			int in  = rndUniform(rndGen) * (nodes.size() - nOutputs);
			int out = rndUniform(rndGen) * (nodes.size() - nInputs);

			if(nodes[in].type == Output) in  += nOutputs;
			if(nodes[out].type == Input) out += nInputs;

            int timeoutTries = nodes.size() * nodes.size();

			// Make sure the connection doesn't already exist
			for(int i=0; i<connections.size(); i++) {
				if(connections[i].in == in && connections[i].out == out) {
                    if (!timeoutTries--) return; // For n nodes, give up after n^2 tries

                    int in = rndUniform(rndGen) * (nodes.size() - nOutputs);
                    int out = rndUniform(rndGen) * (nodes.size() - nInputs);
                    if(nodes[in].type == Output) in += nOutputs;
                    if(nodes[out].type == Input) out += nInputs;
					i = 0;
				}
			}
			float weight = rndNormal(rndGen);

			// Add new connection (with random weight) to list
			connections.push_back( Connection(in, out, weight, true) );
		}
	}

	// Mutate structure by adding node
	void mutateAddNode() {
        if (connections.size() > 0) {
            // Select and disable random connection
            uniform_int_distribution<int> rndInt(0, connections.size() - 1);
            int iRndCon = rndInt(rndGen);
            connections[iRndCon].expressed = false;

            nodes.push_back(Node(Hidden)); // Add new node to list
            int iNewNode = nodes.size() - 1; // Get new node's index

            // Connect new node between previously connected nodes
            connections.push_back(Connection(
                connections[iRndCon].in,
                iNewNode,
                1.0f,
                true
                ));
            connections.push_back(Connection(
                iNewNode,
                connections[iRndCon].out,
                connections[iRndCon].weight,
                true
                ));
        }
	}

	// Propagate values through network one step
	void computeOneCycle(){
		for(int i=0; i<connections.size(); i++) {
			if(connections[i].expressed) {
				nodes[connections[i].out].preVal += 
					nodes[connections[i].in].postVal * 
					connections[i].weight;
			}
		}
		for(int i=0; i<nodes.size(); i++) {
			nodes[i].postVal = nodes[i].activationFunction(nodes[i].preVal);
            nodes[i].preVal = 0.0f; // TODO: why?
		}
	}

public:
	// Default constructor, needed to initialize array
	NerveSystem() {}
	
	// Create nervous system, specifying number of input and output nodes
	NerveSystem(int nIn, int nOut) {
		normal_distribution<float> rndNormal(0.0f, 1.0f);
		nInputs = nIn;
		nOutputs = nOut;
		// Add input nodes
		for(int i=0; i<nIn; i++) {
			nodes.push_back(Node(Input));
		}
		// Add output nodes and connect them to each input node
		for(int i=0; i<nOut; i++) {
			nodes.push_back(Node(Output));
            int out = i;
            for (int j = 0; j<nIn; j++) {
                int in = j;
                float weight = rndNormal(rndGen);
                connections.push_back(Connection(in, out, weight, true));
            }
		}
	}

    void updateInputs(int nIn) {
        normal_distribution<float> rndNormal(0.0f, 1.0f);
        while (nInputs < nIn) {
            nodes.push_back(Node(Input));
            nInputs++;
            int in = nodes.size()-1; //Last index
            for (int n = 0; n < nodes.size()-1; n++) {
                if (nodes[n].type == Output) {
                    float weight = rndNormal(rndGen);
                    connections.push_back(Connection(in, n, weight, true));
                }
            }
        }

    }

	// Input data into the network and get output
	vector<float> getOutput(vector<float> input) {
		// Clear previous values
//		for(int i=0; i<nodes.size(); i++) {
//			nodes[i].preVal = nodes[i].postVal = 0.0f;
//		}

        int iInput = 0;
        for (int iNode = 0; iNode<nodes.size(); iNode++) {
            if (nodes[iNode].type == Input){
                // It is possible for a sensor cell to die, in that case
                // its signal value will be 0.0f
                try {
                    nodes[iNode].postVal = input.at(iInput);
                }
                catch (const std::out_of_range& oor){
                    nodes[iNode].postVal = 0.0f;
                }
                iInput++;
            }
		}
		for(int i=0; i<nActivationCycles; i++) {
			computeOneCycle();
		}
		vector<float> output;
		for(int i=0; i<nodes.size(); i++) {
			if(nodes[i].type == Output){
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
		if(rndNormal(rndGen) < MUTATION_PROB)
			mutateAddConnection();
		if(rndNormal(rndGen) < MUTATION_PROB)
			mutateAddNode();
	}
	
	string toJSON() {
		string vertices = "";
		string links = "";
		string weights = "";

		bool first = true;
		for(int i=0; i<nodes.size(); i++) {
			if(first) first = false;
			else vertices += ",";
			
			vertices += 
				"{\"i\":" 				 +
				to_string(i) 			 +
				",\"type\":" 			 +
				to_string(nodes[i].type) +
				"}";
		}
		
		first = true;
		for(int i=0; i<connections.size(); i++) {
            if (connections[i].expressed) {
                if (first) first = false;
                else {
                    links += ",";
                    weights += ",";
                }

                links +=
                    "{\"i\":" +
                    to_string(connections[i].in) + ",\"o\":" +
                    to_string(connections[i].out) +
                    "}";
                weights +=
                    to_string(connections[i].weight);
            }
		}
		
		return
			string("\"nervesystem\":{") +
			string("\"vertices\":[")    + vertices + string("],") +
			string("\"links\":[")       + links    + string("],") +
			string("\"weights\":[")     + weights  + string("]}");
	}

	void printMathematica() {
		string vertices = "";
		string links = "";
		string vertexStyle = "";
		string edgeWeights = "";
		for(int i=0; i<nodes.size(); i++) {
			vertices += to_string(i) +
				(i==nodes.size()-1 ? "" : ", ");
			switch(nodes[i].type) {
				case Input : vertexStyle += (to_string(i) + "->Blue, "); break;
				case Output: vertexStyle += (to_string(i) + "->Red, ");  break;
			}
		}
		for(int i=0; i<connections.size(); i++) {
			links +=
				to_string(connections[i].in)  + "->" +
				to_string(connections[i].out) +
				(i==connections.size()-1 ? "" : ", ");
			edgeWeights += to_string(connections[i].weight) +
				(i==connections.size()-1 ? "" : ", ");
		}
		cout<<
            "Graph[{"<<vertices<<"},{"<<links<<"}, "<<
			"VertexStyle->{"<<vertexStyle<<"}, "<<
			"EdgeWeight->{"<<edgeWeights<<"}, "<<
			"VertexLabels->\"Name\"]"<<endl;
	}
};

int NerveSystem::currInnovNumber = 0;

#endif