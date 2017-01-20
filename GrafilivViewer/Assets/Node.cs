using System.Collections;
using System.Collections.Generic;
using UnityEngine;

//http://answers.unity3d.com/questions/12322/drag-gameobject-with-mouse.html
public class Node : MonoBehaviour {
    public void OnDrag()
    {
        transform.position = Input.mousePosition;
    } 
 }
