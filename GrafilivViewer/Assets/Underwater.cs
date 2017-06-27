using UnityEngine;
using System.Collections;
 
public class Underwater : MonoBehaviour {
 
	//This script enables underwater effects. Attach to main camera.
 
    //Define variable
    public int underwaterLevel = 100;
 
    //The scene's default fog settings
    private bool defaultFog = RenderSettings.fog;
    private Color defaultFogColor = RenderSettings.fogColor;
    private float defaultFogDensity = RenderSettings.fogDensity;
    private Material defaultSkybox = RenderSettings.skybox;
    private Material noSkybox;
 
    void Start () {
	    //Set the background color
	    //camera.backgroundColor = new Color(0, 0.4f, 0.7f, 1);
    }
 
    void Update () {
        if (transform.position.y < underwaterLevel)
        {
            RenderSettings.fog = true;
            RenderSettings.fogColor = new Color(104.0f/255, 131.0f/255, 139.0f/255, 0.5f);
            RenderSettings.fogDensity = 0.01f;
            RenderSettings.skybox = noSkybox;
        }
        else
        {
            RenderSettings.fog = defaultFog;
            RenderSettings.fogColor = defaultFogColor;
            RenderSettings.fogDensity = defaultFogDensity;
            RenderSettings.skybox = defaultSkybox;
        }
    }
}