#ifndef GENOME_H
#define GENOME_H

#define MUTATION_PROB 0.7f

#define RANDOM (float)rand()/RAND_MAX

using namespace std;

class Genome {
	enum NodeType {Input, Hidden, Output};
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
			type(t), preVal(0.0f), postVal(0.0f),
			f(ActivationFunction((int)(RANDOM * N_ACTIVATION_FUNCTIONS)))
		{}

		float activationFunction(float x) {
			switch(f) {
				case Sine:	return sin(x);
				case Abs : 	return clamp(1.0f - abs(x), -1, 1);
				case Gaus:  return exp(-(x*x)/2);
				case Id  :	return x;
				case Mod :	return fmod(x, 1);
				default  : 
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
	int nActivationCycles = 1;

	// Mutate each connection with small pertubation
	void mutateConnections() {
		for(int i=0; i<connections.size(); i++)
			if(RANDOM < MUTATION_PROB)
				connections[i].weight *= rnd_normal();
	}

	// Mutate structure by adding connection
	void mutateAddConnection() {
		// If there are connection left to add
		if (connections.size() < nodes.size()*nodes.size()) {
			int in  = RANDOM * (nodes.size()-nOutputs);
			int out = RANDOM * (nodes.size()-nInputs);
			
			if(nodes[in].type == Output) in+=nOutputs;
			if(nodes[out].type == Input) out+=nInputs;

			// Make sure the connection doesn't already exist
			for(int i=0; i<connections.size(); i++) {
				if( connections[i].in == in && connections[i].out == out) {
					printf("in=%d, out=%d already exists, trying another connection\n", in, out);
					in  = connections.size() * RANDOM;
					out = connections.size() * RANDOM;
					i = 0;
				}
			}
			float weight = rnd_normal();
			printf("Adding connection: in=%d, out=%d (weight=%f)\n", in, out, weight);
			// Add new connection (with random weight) to list		
			connections.push_back( Connection(in, out, weight, true) );
		}
	}

	// Mutate structure by adding node
	void mutateAddNode() {
		// Select and disable random connection
		int iRndCon = RANDOM * connections.size();
		connections[iRndCon].expressed = false;

		nodes.push_back(Node(Hidden)); // Add new node to list
		int iNewNode = nodes.size()-1; // Get new node's index
		
		printf("Adding new node (#%d), between #%d and #%d, in place of connection #%d\n",
			iNewNode, connections[iRndCon].in, connections[iRndCon].out, iRndCon
		);
		
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
	
	// Propagate values through network one step
	void computeOneCycle(){
		for(Connection c : connections) {
			if(c.expressed) {
				nodes[c.out].preVal += nodes[c.in].postVal * c.weight;
			/*	printf("\nn[%i].preval += n[%i].postVal * weight\n(%.2f += %.2f * %.2f) = %.2f\n",
					c.out, c.in,
					(nodes[c.out].preVal - nodes[c.in].postVal * c.weight),
					nodes[c.in].postVal, c.weight, nodes[c.out].preVal
				);
			*/
			}
		}
		for(int i=0; i<nodes.size(); i++) {
			nodes[i].postVal = nodes[i].activationFunction(nodes[i].preVal);
		/*	if(nodes[i].type == Output) printf("output:");
			printf("( %.2f %.2f )\t", nodes[i].preVal, nodes[i].postVal);
		*/	nodes[i].preVal = 0.0f; // TODO: why?
		}
		//printf("\n");
	}

public:
	// Create genome, specifying number of input and output nodes
	Genome(int nIn, int nOut) {
		srand(time(0));
		
		nInputs = nIn;
		nOutputs = nOut;
		// Add input nodes
		for(int i=0; i<nIn; i++) {
			nodes.push_back(Node(Input));
		}
		// Add output nodes and connect them each to a random input node
		for(int i=0; i<nOut; i++) {
			nodes.push_back(Node(Output));
			
			int in  = RANDOM * nIn;
			int out = nodes.size()-1;
			float weight = rnd_normal();
		/*	printf("\nConnecting node %d to node %d (weight=%.2f)\n",
				in, out, weight
			);
		*/	connections.push_back(Connection(in, out, weight, true));
		}
	}
	
	// Input data into the network and get output
	vector<float> getOutput(vector<float> input) {
		// Clear previous values
		for(Node n : nodes) {
			n.preVal = n.postVal = 0.0f;
		}
		
		// Set input (assuming they are located first in the array)
		for(int i=0; i<input.size(); i++) {
			if(nodes[i].type == Input)
				nodes[i].postVal = input[i];
			else {
				printf( "Node %d is not input, has type %i", i, nodes[i].type);
				printGenome();
			}
		}

		for(int i=0; i<nActivationCycles; i++) {
			computeOneCycle();
		}

		vector<float> output;

		for(Node n : nodes) {
			if(n.type == Output){
				//TODO: should we use activationFunction here?
				output.push_back(n.activationFunction(n.postVal));
			}
		}

		return output;
	}
	
	// Mutate network
	void mutate(){
		mutateConnections();
		if(RANDOM < MUTATION_PROB)
			mutateAddConnection();
		if(RANDOM < MUTATION_PROB)
			mutateAddNode();
	}
	
	// Print genome in human readable format
	void printGenome() {
		printf("\n\nGenome: \n");
		printf("%d links:\t", connections.size());
		for(Connection c : connections) {
			printf("#%d:[%d-(w%.2f)%c>%d] ",
				c.innovNumber, c.in, c.weight, c.expressed ? '-' : '/', c.out
			);
		}
		printf("\n%d nodes:\t", nodes.size());
		for(int i=0; i<nodes.size(); i++) {
			char type;
			switch(nodes[i].type) { 
				case Input : type='i'; break;
				case Output: type='o'; break;
				case Hidden: type='h'; break;
			}
			printf("%d:[%c, f%i]\t", i, type, nodes[i].f);
		}
	}
	void printMathematica() {
		string vertices = "";
		string links = "";
		string vertexStyle = "";
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
		}
		cout<< "Graph[{"<<vertices<<"},{"<<links<<"}, VertexStyle->{"<<vertexStyle<<"}, VertexLabels->\"Name\"]"<<endl;
	}
};

int Genome::currInnovNumber = 0;

#endif