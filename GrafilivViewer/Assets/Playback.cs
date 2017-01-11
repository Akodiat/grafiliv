using System.Collections;
using System.Collections.Generic;
using UnityEngine.UI;
using UnityEngine;

public class Playback : MonoBehaviour {
    [SerializeField]
    private InputField frameField;

    [SerializeField]
    private Toggle toggle;


	// Use this for initialization
	void Start () {
		
	}
	
	// Update is called once per frame
	void Update () {
        if (toggle.isOn)
        {
            int frame = int.Parse(frameField.text);
            frame++;
            frameField.text = ""+frame;
        }
	}
}
