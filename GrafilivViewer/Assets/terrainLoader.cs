using System.Collections;
using System.Collections.Generic;
using System.IO;
using System.Text.RegularExpressions;
using UnityEngine;

public class terrainLoader : MonoBehaviour
{

    // Use this for initialization
    void Start()
    {
        string line;

        int nEggs = 0;
        StreamReader sr = new StreamReader("conf.txt");
        
        while (sr.Peek() >= 0) 
        { 
            try
            {
                line = sr.ReadLine();
            }
            catch (FileNotFoundException)
            {
                return;
            }           
        
            //Disregard everything after the semicolon:
            string s = line.Split(';')[0];

            //Remove spaces:
            s = Regex.Replace(s, @"\s+", "");

            //Disregard empty lines:
            if (s.Length == 0) continue;

            string[] pair = s.Split('=');
            string key = pair[0];
            string val = pair[1];

            if (key.Equals("w"))
            {
                string[] ws = val.Split(',');
                Vector3 w = new Vector3(
                    int.Parse(ws[0]),
                    int.Parse(ws[1]), 
                    int.Parse(ws[2])
                );
                print(w);
                Vector3 oldSize = this.transform.localScale;
                this.transform.localScale *= Mathf.Min(w.x, w.z) / 9.0f;
                Vector3 newSize = this.transform.localScale;

                Vector3 posDiff = newSize - oldSize;
                //posDiff.z += oldSize.z;

                this.transform.position = new Vector3(w.x / 2 - posDiff.x, -posDiff.y, w.z / 2 + posDiff.z/2);

                break;
            }
        }
    }

    // Update is called once per frame
    void Update()
    {

    }
}
