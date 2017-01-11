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

    public GameObject photocyte;
    public GameObject phagocyte;

    [SerializeField]
    private InputField frameField;

    private int frame;

    private List<GameObject> particles = new List<GameObject>();

    private FileInfo[] files;

	// Use this for initialization
	void Start () {
        DirectoryInfo dir = new DirectoryInfo("output");
        files = dir.GetFiles();
    }
    // Update is called once per frame
    void Update () {
        int newFrame = int.Parse(frameField.text);

        if (frame != newFrame)
        {
            frame = newFrame;
            particles.ForEach(particle => Destroy(particle));

            string line;
            
            try {
                line = files[frame].OpenText().ReadToEnd();
                /*
                using (StreamReader sr = new StreamReader("output/frame" + frame + ".json")) {
                    line = sr.ReadToEnd();
                }
                */
            }
            catch (FileNotFoundException e) {
                return;
            }


            Particle[] p = JsonHelper.FromJson<Particle>(line);

            //print("Number of particles: " + p.Length);

            for (int i = 0; i < p.Length; i++)
            {
                if (p[i].pt == ParticleType.Cell)
                {
                    Vector3 position = new Vector3(
                        p[i].x, p[i].y, p[i].z
                    );
                    switch (p[i].ct)
                    {
                        case CellType.Photo:
                            particles.Add(Instantiate(photocyte, position, Quaternion.identity));
                            break;
                        case CellType.Digest:
                            particles.Add(Instantiate(phagocyte, position, Quaternion.identity));
                            break;
                    }
                }
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