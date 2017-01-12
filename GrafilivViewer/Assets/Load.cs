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
    public Material mPhoton;

    [SerializeField]
    private InputField frameField;

    private int frame;

    private List<GameObject> particles = new List<GameObject>();

    private FileInfo[] files;

	// Use this for initialization
	void Start () {
        DirectoryInfo dir = new DirectoryInfo("output");
        files = dir.GetFiles();
        frame = 0;

        string line;   
        try {
            line = files[frame].OpenText().ReadToEnd();
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
                switch (p[i].ct)
                {
                    case CellType.Photo:
                        renderer.material = mPhotocyte;
                        break;
                    case CellType.Digest:
                        renderer.material = mPhagocyte;
                        break;
                }
            }
            else if (p[i].pt == ParticleType.Energy)
                renderer.material = mPhoton;
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
                line = files[frame].OpenText().ReadToEnd();
            }
            catch (ArgumentOutOfRangeException)
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
                    switch (p[i].ct)
                    {
                        case CellType.Photo:
                            renderer.material = mPhotocyte;
                            break;
                        case CellType.Digest:
                            renderer.material = mPhagocyte;
                            break;
                    }
                }
                else if (p[i].pt == ParticleType.Energy)
                    renderer.material = mPhoton;
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