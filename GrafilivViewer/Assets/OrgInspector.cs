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
    public Texture2D point;

    private List<Node> nodes = new List<Node>();
    private List<Connection> connections = new List<Connection>();

    private static int trimRadius = 25;
    private static int springDistance = 300;

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
        //this.loadOrganism(0);
	}
	
	// Update is called once per frame
	void Update () {
	}

    //Remove old nodes and connections
    public void unloadOrganism()
    {
        foreach (Node node in nodes)
            Destroy(node.getGameObject());
        nodes.Clear();
        connections.Clear();
    }

    void loadOrganism(int organismID)
    {
        unloadOrganism();

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

        Organism org = UnityEngine.JsonUtility.FromJson<Organism>(line);

        for (int i = 0; i < org.genome.vertices.Count; i++)
        {
            var v = org.genome.vertices[i];
            Node node = new Node(i, v.type, v.f, nodePrefab, canvas);
            nodes.Add(node);
        }
        for(int i=0; i < org.genome.links.Count; i++)
        {
            var c = org.genome.links[i];
            if (c.i == c.o) continue; //Ignore self connections at the moment
            float w = org.genome.weights[i];
            connections.Add(new Connection(nodes[c.i], nodes[c.o], w));
        }
    }

    void OnGUI() {
        if (Event.current == null)
            return;

        if (Event.current.type != EventType.repaint)
            return;
        foreach (Connection c in connections)
        {
            Vector2 p1 = c.getNode1().getPosition();
            Vector2 p2 = c.getNode2().getPosition();
            Vector2 dx = p2 - p1;
            DrawLine(
                p1 + (dx.normalized * trimRadius),
                p2 - (dx.normalized * trimRadius),
                (int) (10 * c.getWeight())
            );
        }
    }

    private void DrawLine(Vector2 start, Vector2 end, int width)
    {
        Vector2 d = end - start;
        float a = Mathf.Rad2Deg * Mathf.Atan(d.y / d.x);
        if (d.x < 0)
            a += 180;

        int width2 = (int)Mathf.Ceil(width / 2);

        GUIUtility.RotateAroundPivot(a, start);
        GUI.DrawTexture(new Rect(start.x, start.y - width2, d.magnitude, width), point);
        GUIUtility.RotateAroundPivot(-a, start);
    }

    private class Node
    {
        private GameObject gameObject;
        private NodeType type;
        private ActivationFunction activationFunction;
        private int id;

        public Node(int id, NodeType type, ActivationFunction f, GameObject nodePrefab, GameObject canvas)
        {
            GameObject o = Instantiate(nodePrefab) as GameObject;
            o.transform.SetParent(canvas.transform);
            Color color;

            float x = 500 + id*10;
            float y;
            switch (type)
            {
                case NodeType.Input:
                    color = Color.blue;
                    y = 800;
                    break;
                case NodeType.Hidden:
                    color = Color.gray;
                    y = 500;
                    break;
                case NodeType.Output:
                    color = Color.green;
                    y = 200;
                    break;
                default: color = Color.black;
                    y = 500;
                    break;
            }
            o.GetComponent<Image>().color = color;
            o.transform.position = new Vector3(x, y);

            this.id = id;
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
    public int parent;
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

