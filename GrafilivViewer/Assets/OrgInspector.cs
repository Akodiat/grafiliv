using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UI;
using System.IO;

public enum NodeType {Input, Hidden, Output};
public enum ActivationFunction {
	Sine, Abs, Id, Mod, Gaus
};

public class OrgInspector : MonoBehaviour {

    public GameObject nodePrefab;
    public GameObject canvas;

    private List<Node> nodes = new List<Node>();
    private List<Connection> connections = new List<Connection>();

    private static int trimRadius = 25;
    private static int springDistance = 100;

	// Use this for initialization
	void Start () {
        /*
        Node node1 = new Node(Node.NodeType.Input, Node.ActivationFunction.Abs, nodePrefab, canvas);
        Node node2 = new Node(Node.NodeType.Output, Node.ActivationFunction.Abs, nodePrefab, canvas);
        Node node3 = new Node(Node.NodeType.Output, Node.ActivationFunction.Abs, nodePrefab, canvas);

        nodes.Add(node1); nodes.Add(node2); nodes.Add(node2);

        connections.Add(new Connection(node1, node2, 3.0f));
        connections.Add(new Connection(node1, node3, 2.0f));
        */
        this.loadOrganism(0);
	}
	
	// Update is called once per frame
	void Update () {
	}

    void loadOrganism(int organismID)
    {
        nodes.Clear();
        connections.Clear();

        string line;
        try
        {
            StreamReader sr = new StreamReader("organisms/org" + organismID + ".json");
            line = sr.ReadToEnd();
        }
        catch (FileNotFoundException)
        {
            return;
        }

        print("File loaded: " + line);

        Organism org = UnityEngine.JsonUtility.FromJson<Organism>(line);
        print("Radius x: " + org.genome.radius[0]);

        print("Weight of link 0: " + org.genome.weights[0]);

        print(org.genome.links.Count + " links");
        print(org.nervesystem.links.Count + " nerveLinks");
    }

    void OnGUI() {
        foreach (Connection c in connections)
        {
            Vector2 p1 = c.getNode1().getPosition();
            Vector2 p2 = c.getNode2().getPosition();
            Vector2 dx = p2 - p1;
            Drawing.DrawLine(
                p1, //+ (dx.normalized * trimRadius),
                p2, //- (dx.normalized * trimRadius),
                Color.black,
                c.getWeight()
            );
        }
    }

    private class Node
    {
        private GameObject gameObject;
        private NodeType type;
        private ActivationFunction activationFunction;

        public Node(NodeType type, ActivationFunction f, GameObject nodePrefab, GameObject canvas)
        {
            GameObject o = Instantiate(nodePrefab) as GameObject;
            o.transform.SetParent(canvas.transform);
            o.transform.position = new Vector3(500 + Random.value, 500 + Random.value);
            
            this.gameObject = o;
            this.type = type;
            this.activationFunction = f;
        }

        public Vector2 getPosition()
        {
            Vector2 pos = gameObject.GetComponent<Rigidbody2D>().position;
            pos.y = Screen.height - pos.y;
            return pos;
        }

        public GameObject getGameObject() { return gameObject; }
        public NodeType getType() { return type; }
        public ActivationFunction getActivationFunction() { return activationFunction; }
    }

    private class Connection
    {
        private Node node1, node2;
        private float weight;

        public Connection(Node node1, Node node2, float weight)
        {
            this.node1 = node1;
            this.node2 = node2;
            this.weight = weight;

            SpringJoint2D spring = node1.getGameObject().AddComponent<SpringJoint2D>();
            spring.connectedBody = node2.getGameObject().GetComponent<Rigidbody2D>();
            spring.distance = springDistance;
        }

        public Node getNode1() { return node1; }
        public Node getNode2() { return node2; }
        public float getWeight() { return weight; }
    }
}

[System.Serializable]
public class Organism
{
    public Genome genome;
    public NerveSystem nervesystem;
}

[System.Serializable]
public struct Genome
{
    public List<GenomeVert> vertices;
    public List<Link> links;
    public List<float> weights;
    public int[] radius;
}

[System.Serializable]
public struct NerveSystem
{
    public List<NerveVert> vertices;
    public List<Link> links;
    public List<float> weights;
}

[System.Serializable]
public struct Link
{
    public int i, o;
}

[System.Serializable]
public struct GenomeVert
{
    public int i;
    public NodeType type;
    public ActivationFunction f;
}

[System.Serializable]
public struct NerveVert
{
    public int i;
    public NodeType type;
}

