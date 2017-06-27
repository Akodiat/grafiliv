using UnityEngine;
using System.Collections;

// https://forum.unity3d.com/threads/replacing-spheres-by-sprites-if-they-are-out-camera-range-3d-space-flier.299743/
// Makes object (only tested for spheres) too close to the Camera.farClipPlane invisible
// and replaces it with a sprite at the edge of the farClipPlane.

// if anyone finds this code usefull just take it :)

public class SeeClose : MonoBehaviour
{

    public GameObject FarObject; // the sprite to be used for replacing the object if too far

    private GameObject farObject; //Local instance of the sprite
    private Camera myCamera; // the scenes camera

    // Use this for initialization
    void Start()
    {

        //Find the camera
        myCamera = Camera.FindObjectOfType<Camera>();

        if (!myCamera)
        {
            Debug.Log("SeeClose needs a Camera component!");
            enabled = false;
            return;
        }

        if (!FarObject)
        {
            Debug.Log("SeeClose needs a FarObject set for " + name);
            enabled = false;
            return;
        }

        //Instantiate a sprite for replacement
        farObject = Instantiate(FarObject) as GameObject;
        farObject.GetComponent<Renderer>().enabled = false;
        farObject.name = this.name + "_farObject";
    }


    // Update is called once per frame
    void Update()
    {

        if (isClose())
        {
            // when object is close render the object
            GetComponent<Renderer>().enabled = true;
            farObject.GetComponent<Renderer>().enabled = false;
        }
        else
        {
            // when object far render the sprite
            GetComponent<Renderer>().enabled = false;

            // get the direction of the object
            Vector3 heading = transform.position - myCamera.transform.position;
            float distance = heading.magnitude;
            Vector3 direction = heading / distance; // This is now the normalized direction.

            // Find the position of the sprite in camera range
            float fakeDistance = 19000;
            Vector3 pos = direction * fakeDistance;
            // find the scale the sprite must be
            Vector3 scale = transform.localScale * (fakeDistance / distance);

            // Remove even the sprites if object is too small
            if ((scale.x < 100) && (scale.y < 100) && (scale.z < 100))
            {
                farObject.GetComponent<Renderer>().enabled = false;
                return;
            }

            // make the sprite always look at the player
            farObject.transform.LookAt(myCamera.transform);

            // Set position and scale for the sprite
            farObject.transform.position = pos;
            farObject.transform.localScale = scale;

            // Sort the sprite according to the distance to the player
            farObject.GetComponent<Renderer>().sortingOrder = int.MaxValue - (int)(distance / 100);
            farObject.GetComponent<Renderer>().enabled = true;
        }
    }

    private bool isClose()
    {
        // see if the distance of the sphere is inside the farClipPlane

        // Since some objects can be very large I need to remove it's scale from
        // the farClipPlane

        // I use sqrMagnitude for speed.

        float CameraViewDist = (myCamera.farClipPlane - transform.localScale.x) / 2; //I don't understand why i need /2 here but it seems to work.
        float FarDist = (myCamera.transform.position - transform.position).sqrMagnitude;

        return CameraViewDist * CameraViewDist > FarDist;
    }
}