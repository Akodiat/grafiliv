using System;
using UnityEngine;
using System.IO;
using UnityEngine.UI;
using System.Collections.Generic;

public enum CellType
{
    Photo, Digest, Sting, Vascular, Fat, Sense, Egg,
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


    [SerializeField]
    private InputField frameField;

    private int frame;

    private List<GameObject> particles = new List<GameObject>();

    private FileInfo[] files;
    private DirectoryInfo dir;

    
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
            default:
                return null;
        }
    }

	// Use this for initialization
	void Start () {
        //dir = new DirectoryInfo("output");
        //files = dir.GetFiles();
        frame = 0;

        string line;   
        try {
            StreamReader sr = new StreamReader("output/frame" + frame + ".json");
            line = sr.ReadToEnd();
        }
        catch (FileNotFoundException) {
            return;
        }
        Particle[] p = JsonHelper.FromJson<Particle>(line);

        for (int i = 0; i < p.Length; i++)
        {

            Vector3 position = new Vector3(
                p[i].x, p[i].y, p[i].z
            );

            GameObject particle = Instantiate(particlePrefab, position, Quaternion.identity);
            MeshRenderer renderer = particle.GetComponent(typeof(MeshRenderer)) as MeshRenderer;
                
            if (p[i].pt == ParticleType.Cell)
            {
                renderer.material = getMaterial(p[i].ct);
            }
            particles.Add(particle);
        }
    }
    // Update is called once per frame
    void Update () {
        int newFrame = int.Parse(frameField.text);

        if (frame != newFrame)
        {
            frame = newFrame;

            string line;

            try
            {
                StreamReader sr = new StreamReader("output/frame" + frame + ".json");
                line = sr.ReadToEnd();
            }
            catch (FileNotFoundException)
            {
                return;
            }

            Particle[] p = JsonHelper.FromJson<Particle>(line);

            int i = 0;
            while (i < p.Length)
            {
                Vector3 position = new Vector3(
                        p[i].x, p[i].y, p[i].z
                );

                //print(i);

                if (i >= particles.Count)
                {
                    GameObject particle = Instantiate(particlePrefab, position, Quaternion.identity);
                    particles.Add(particle);
                }
                else
                {
                    particles[i].transform.position = position;
                }

                MeshRenderer renderer = particles[i].GetComponent(typeof(MeshRenderer)) as MeshRenderer;
                if (p[i].pt == ParticleType.Cell)
                {
                    renderer.material = getMaterial(p[i].ct);
                }

                i++;
            }
            while (i < particles.Count)
            {
                Destroy(particles[i]);
                particles.RemoveAt(i);
            }
        }
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
    /*
    public string playerId;
    public string playerLoc;
    public string playerNick;
    */
    public ParticleType pt;
    public CellType ct;
    public int o;
    public float x;
    public float y;
    public float z;
}