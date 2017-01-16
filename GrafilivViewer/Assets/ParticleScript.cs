using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UI;

public class ParticleScript : MonoBehaviour {

    public GameObject organismIndicator;

    private static GameObject infoBox;

    private int organismId;

    void setOrgId(int o)
    {
        organismId = o;
    }

	// Use this for initialization
	void Start () {
        infoBox = GameObject.Find("Content");
	}

    void OnMouseDown()
    {
        //print("You clicked on me!!!");
        Text text = infoBox.GetComponent("Text") as Text;
        foreach (GameObject g in GameObject.FindGameObjectsWithTag("SelectedOrganism"))
        {
            Destroy(g);
        }
        Instantiate(organismIndicator, this.transform.position, Quaternion.identity);
        //print(text.text);
        text.text = "Organism ID: " + organismId;
    }
}
