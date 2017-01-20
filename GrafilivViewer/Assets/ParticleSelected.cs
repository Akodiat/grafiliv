using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UI;

public class ParticleSelected : MonoBehaviour {

    public GameObject organismIndicator;

    private static GameObject infoBox;
    private static GameObject organismInspector;

    private Particle p;

    void setParticle(Particle p)
    {
        this.p = p;
    }

    // Use this for initialization
    void Start()
    {
        infoBox = GameObject.Find("Infobox");
        organismInspector = GameObject.Find("OrganismInspector");
    }

    void OnMouseDown()
    {
        //print("You clicked on me!!!");
        Text text = infoBox.GetComponent("Text") as Text;
        foreach (GameObject g in GameObject.FindGameObjectsWithTag("SelectedOrganism"))
            Destroy(g);
        Instantiate(organismIndicator, this.transform.position, Quaternion.identity);
        //print(text.text);
        text.text =
            "Particle type: " + p.pt +
            "\nCell type: " + p.ct +
            "\nOrganism ID: " + p.o +
            "\nEnergy: " + p.e;
        organismInspector.SendMessage("loadOrganism", p.o);
    }
}
