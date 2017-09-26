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
                posDiff.z += oldSize.z;

                //this.transform.position = new Vector3(w.x / 2 - posDiff.x, -posDiff.y, w.z / 2 + posDiff.z/2);
                this.transform.position = new Vector3(w.x / 2, -2, w.z / 2);

                GameObject water = GameObject.Find("water");
                water.transform.position = new Vector3(w.x / 2, w.y / 2, w.z / 2);

                break;
            }
        }
        Mesh mesh = GetComponent<MeshFilter>().mesh;
        List<Vector3> oldVertices = new List<Vector3>(mesh.vertices);
        List<Vector3> newVertices = new List<Vector3>();

        sr = new StreamReader("terrain/mesh.csv");      
        line = sr.ReadLine(); //Read heading line so that we don't try to parse it 
        int i = 0;
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
            TerrainVertex p = new TerrainVertex(line);

            newVertices.Add(p.pos);
        }
        print("Old vertices length:");
        print(oldVertices.Count);
        print("New vertices length:");
        print(newVertices.Count);
        //mesh.vertices = vertexList;
    }

    // Update is called once per frame
    void Update()
    {

    }
}

[System.Serializable]
public class TerrainVertex
{
    public TerrainVertex(string s)
    {
        string[] fields = s.Split(',');
        pos = new Vector3(
            float.Parse(fields[0]),
            float.Parse(fields[1]),
            float.Parse(fields[2])
        );
    }
    public Vector3 pos;
}
