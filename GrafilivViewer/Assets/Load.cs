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
    public Material mVascular;


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
            case CellType.Vascular:
                return mVascular;
            default:
                return null;
        }
    }

	// Use this for initialization
	void Start () {
        //dir = new DirectoryInfo("output");
        //files = dir.GetFiles();
        frame = 0;

        StreamReader sr;

        string line;   
        try {
            sr = new StreamReader("output/frame" + frame + ".json");
        }
        catch (FileNotFoundException) {
            return;
        }
        while (sr.Peek() >= 0)
        {
            line = sr.ReadLine();

            Particle p = Particle.parseParticle(line);

            Vector3 position = new Vector3(
                p.x, p.y, p.z
            );
            GameObject particle = Instantiate(particlePrefab, position, Quaternion.identity);
            particle.SendMessage("setOrgId", p.o);

            MeshRenderer renderer = particle.GetComponent(typeof(MeshRenderer)) as MeshRenderer;

            if (p.pt == ParticleType.Cell)
            {
                renderer.material = getMaterial(p.ct);
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

            StreamReader sr;

            string line;
            try
            {
                sr = new StreamReader("output/frame" + frame + ".json");
            }
            catch (FileNotFoundException)
            {
                return;
            }

            int i = 0;
            while (sr.Peek() >= 0)
            {
                line = sr.ReadLine();

                Particle p = Particle.parseParticle(line);

                Vector3 position = new Vector3(p.x, p.y, p.z);

                //print(i);
                /*
                if (i >= particles.Count)
                {
                    GameObject particle = Instantiate(particlePrefab, position, Quaternion.identity);
                    particles.Add(particle);
                }
                else
                {*/
                    particles[i].transform.position = position;
                //}

                particles[i].SendMessage("setOrgId", p.o);

                MeshRenderer renderer = particles[i].GetComponent(typeof(MeshRenderer)) as MeshRenderer;
                if (p.pt == ParticleType.Cell)
                {
                    renderer.material = getMaterial(p.ct);
                }

                i++;
            }
            /*
            while (i < particles.Count)
            {
                Destroy(particles[i]);
                particles.RemoveAt(i);
            }
            */
        }
    }
}

public class Particle
{
    public ParticleType pt;
    public CellType ct;
    public int o;
    public float x;
    public float y;
    public float z;

    public Particle(
        ParticleType pt,
        CellType ct,
        int o,
        float x,
        float y,
        float z
    )
    {
        this.pt = pt; this.ct = ct; this.o = o;
        this.x = x; this.y = y; this.z = z;
    }

    public static Particle parseParticle(string line)
    {
        string[] ps = line.Split(',');
        Particle p = new Particle(
            (ParticleType)int.Parse(ps[0]),
            (CellType)int.Parse(ps[1]),
            int.Parse(ps[2]),
            float.Parse(ps[3]),
            float.Parse(ps[4]),
            float.Parse(ps[5])
            );
        return p;
    }
}