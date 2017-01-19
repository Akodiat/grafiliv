using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UI;

public class OrgInspector : MonoBehaviour {

    public GameObject nodePrefab;
    public GameObject canvas;

    private Node node1;
    private Node node2;
    private Connection connection;

	// Use this for initialization
	void Start () {
        node1 = new Node(Node.NodeType.Input, Node.ActivationFunction.Abs, nodePrefab, canvas);
        node2 = new Node(Node.NodeType.Output, Node.ActivationFunction.Abs, nodePrefab, canvas);

        connection = new Connection(node1, node2, 3.0f);
	}
	
	// Update is called once per frame
	void Update () {
	}

    void OnGUI() {
        Drawing.DrawLine(
            connection.getNode1().getPosition(),
            connection.getNode2().getPosition(),
            Color.black,
            connection.getWeight()
        );
    }

    private class Node
    {
        public enum NodeType {Input, Hidden, Output};
        public enum ActivationFunction {
		    Sine, Abs, Id, Mod, Gaus
	    };

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
            spring.distance = 100;
        }

        public Node getNode1() { return node1; }
        public Node getNode2() { return node2; }
        public float getWeight() { return weight; }
    }
}