using System.Collections;
using System.Collections.Generic;
using UnityEngine.UI;
using UnityEngine;

public class Playback : MonoBehaviour {
    [SerializeField]
    private InputField frameInput;

    [SerializeField]
    private InputField frameSkipInput;

    [SerializeField]
    private Toggle toggle;


	// Use this for initialization
	void Start () {
		
	}
	
	// Update is called once per frame
	void Update () {
        if (toggle.isOn)
        {
            nextFrame();
        }
	}

    public void nextFrame()
    {
        int frame = int.Parse(frameInput.text);
        int skip = int.Parse(frameSkipInput.text);
        frame += skip;
        frameInput.text = "" + frame;
    }
}
