using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UI;

public class ParticleSelected : MonoBehaviour {

    public GameObject particleIndicator;
    public GameObject organismIndicator;

    private static GameObject infoBox;
    private static GameObject organismInspector;

    private Particle p;
    private List<Vector3> orgPositions;

    void setParticle(Particle p)
    {
        this.p = p;
    }

    void setOrganism(List<Vector3> orgPositions)
    {
        this.orgPositions = orgPositions;
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
        Instantiate(particleIndicator, this.transform.position, Quaternion.identity);
        foreach (Vector3 pos in orgPositions)
        {
            Instantiate(organismIndicator, pos, Quaternion.identity);
        }

        text.text =
            "Particle type: " + p.pt +
            "\nCell type: " + p.ct +
            "\nOrganism ID: " + p.o +
            "\nEnergy: " + p.e;
        organismInspector.SendMessage("loadOrganism", p.o);
    }
}
