using System;
using UnityEngine;
using System.IO;
using UnityEngine.UI;
using System.Collections.Generic;

public enum CellType
{
    Photo, Digest, Sting, Vascular, Fat, Sense, Egg, Buoyancy,
    N_CELL_TYPES
};
public enum ParticleType
{
    Cell, Energy, Pellet, Buffer,
    N_PARTICLE_TYPES
};

public class Load : MonoBehaviour {

    public GameObject particlePrefab;
    public Material mPhotocyte;
    public Material mPhagocyte;
    public Material mDevorocyte;
    public Material mLipocyte;
    public Material mSensor;
    public Material mEgg;
    public Material mVascular;
    public Material mBouyancy;

    public Material mPellet;

    public Toggle inspectToggle;
    public Toggle recordToggle;

    public Text eggCounter;

    public GameObject terrain;


    [SerializeField]
    private InputField frameField;

    private int frame;

    private List<GameObject> particles = new List<GameObject>();

    private FileInfo[] files;
    private DirectoryInfo recordDir;
    
    private List<Vector3> GetOrgPositions(int orgID)
    {
        /*
        List<Vector3> positions = new List<Vector3>();
        foreach (GameObject p in particles)
            if (p.o == orgID)
                positions.Add(new Vector3(p.x, p.y, p.z));
        return positions;
        */
        return new List<Vector3>();
    }
    
    
    private Material getMaterial(CellType ct){
        switch (ct)
        {
            case CellType.Photo:
                return mPhotocyte;
            case CellType.Digest:
                return mPhagocyte;
            case CellType.Fat:
                return mLipocyte;
            case CellType.Sense:
                return mSensor;
            case CellType.Sting:
                return mDevorocyte;
            case CellType.Egg:
                return mEgg;
            case CellType.Vascular:
                return mVascular;
            case CellType.Buoyancy:
                return mBouyancy;
            default:
                return null;
        }
    }

	// Use this for initialization
	void Start () {
        //dir = new DirectoryInfo("output");
        //files = dir.GetFiles();
        recordDir = Directory.CreateDirectory(Application.dataPath + "/record");

        frame = 0;
    }
    // Update is called once per frame
    void Update () {
        int newFrame = int.Parse(frameField.text);

        if (frame != newFrame)
        {
            frame = newFrame;
            Refresh();
            if (recordToggle.isOn)
                Application.CaptureScreenshot(recordDir.FullName+"/frame"+frame+".png");
        }
    }

    public void Refresh()
    {
        string line;

        int nEggs = 0;
        StreamReader sr = new StreamReader("output/frame" + frame + ".csv");
        
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
            Particle p = new Particle(line);

            Vector3 position = new Vector3(
                    p.x, p.y, p.z
            );

            if (i >= particles.Count)
            {
                GameObject particle = Instantiate(particlePrefab, position, Quaternion.identity);
                particles.Add(particle);
            }
            else
            {
                particles[i].transform.position = position;
            }

            //Only enable collision if inspect is on:
            particles[i].GetComponent<SphereCollider>().enabled = inspectToggle.isOn;
            if (inspectToggle.isOn)
            {
                particles[i].SendMessage("setParticle", p);
                particles[i].SendMessage("setOrganism", GetOrgPositions(p.o));
            }

            MeshRenderer renderer = particles[i].GetComponent(typeof(MeshRenderer)) as MeshRenderer;
            Transform transform = particles[i].GetComponent(typeof(Transform)) as Transform;
            transform.localScale = 2 * new Vector3(p.r, p.r, p.r);

            if (p.pt == ParticleType.Cell)
            {
                renderer.material = getMaterial(p.ct);
                if (p.ct == CellType.Egg)
                    nEggs++;
            }
            else if (p.pt == ParticleType.Pellet)
            {
                renderer.material = mPellet;
            }

            i++;
        }
        while (i < particles.Count)
        {
            Destroy(particles[i]);
            particles.RemoveAt(i);
        }
        eggCounter.text = "# of egg cells: " + nEggs;
    }
}


public static class JsonHelper
{
    public static T[] FromJson<T>(string json)
    {
        Wrapper<T> wrapper = UnityEngine.JsonUtility.FromJson<Wrapper<T>>(json);
        return wrapper.Items;
    }

    public static string ToJson<T>(T[] array)
    {
        Wrapper<T> wrapper = new Wrapper<T>();
        wrapper.Items = array;
        return UnityEngine.JsonUtility.ToJson(wrapper);
    }

    [Serializable]
    private class Wrapper<T>
    {
        public T[] Items;
    }
}

[System.Serializable]
public class Particle
{
    public Particle(string s)
    {
        string[] fields = s.Split(',');
        pt = (ParticleType) int.Parse(fields[0]);
        ct = (CellType)     int.Parse(fields[1]);
        o =                 int.Parse(fields[2]);
        e =               float.Parse(fields[3]);
        r =               float.Parse(fields[4]);
        x =               float.Parse(fields[5]);
        y =               float.Parse(fields[6]);
        z =               float.Parse(fields[7]);
    }
    public ParticleType pt;
    public CellType ct;
    public int o;
    public float e;
    public float r;
    public float x;
    public float y;
    public float z;
}
