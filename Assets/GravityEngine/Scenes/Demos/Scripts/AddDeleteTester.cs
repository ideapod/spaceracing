using UnityEngine;
using System.Collections;
using System.Collections.Generic;

/// <summary>
/// Add/delete tester.
/// Create and remove:
/// - massive, massless, Kepler orbiting bodies (on screen buttons)
/// - binary planets (key B)
/// - add a moon to the first massive planet (key M)
/// 
/// Used to demonstrate dynamic addition and as a test case for dynamic additions during development. 
/// 
/// Objects with trail renders require some special care. The trail renders in the prefabs are disabled and
/// are only enabled once the object has been added to GE AND enough time has gone by for the GE to have set the
/// tranform. 
/// 
/// </summary>
public class AddDeleteTester : MonoBehaviour {

	public GameObject orbitingPrefab;
    public GameObject binaryPrefab;
    public GameObject hyperPrefab;
    public GameObject fixedObject;

    public GameObject star;

	public  float maxRadius= 30f;
	public  float minRadius = 5f;
    public float moonRadius = 2f;
    public float fixedStarRadius = 20f;

    public float maxEccentricity = 0f;
    public float minEccentricity = 0f;
    public float mass = 0.001f;

    //! test flag to do random add/deletes
    public bool runChaosMonkey = false;

	private List<GameObject> massiveObjects; 
	private List<GameObject> masslessObjects; 
	private List<GameObject> keplerObjects;
    private List<GameObject> binaryObjects;
    private List<GameObject> moonObjects;
    private List<GameObject> hyperObjects;
    private List<GameObject> fixedObjects;

    private Color[] colors = { Color.red, Color.white, Color.blue, Color.cyan, Color.gray, Color.green, Color.magenta, Color.yellow};
	private int colorIndex = 0; 

	// Use this for initialization
	void Awake () {
		massiveObjects = new List<GameObject>();
		masslessObjects = new List<GameObject>();
		keplerObjects = new List<GameObject>();
        binaryObjects = new List<GameObject>();
        moonObjects = new List<GameObject>();
        hyperObjects = new List<GameObject>();
        fixedObjects = new List<GameObject>();
    }

    // strings passed in by button functions in UI
    private const string MASSIVE = "massive";
    private const string MASSLESS = "massless";
	private const string KEPLER = "kepler";
    private const string BINARY = "binary";
    private const string MOON = "moon";
    private const string HYPER = "hyper";
    private const string FIXED = "fixed";

    private string[] adTypes;

    void Start() {
        adTypes = new string[] { MASSIVE, MASSLESS, KEPLER, BINARY, MOON, HYPER, FIXED };
    }

    public void AddBody(string bodyType) {
        if (bodyType != MOON) {
            AddBodyToParent(bodyType, star);
        } else {
            // Add a moon to massive body
            if (massiveObjects.Count > 0) {
                AddBodyToParent(bodyType, massiveObjects[massiveObjects.Count - 1]);
            } else {
                Debug.LogWarning("Cannot add moon - no massive objects");
            }
        }
    }


    private void EllipseInit(GameObject go, GameObject parent, string bodyType) {

        OrbitEllipse eb = go.GetComponent<OrbitEllipse>();
        if (eb == null) {
            Debug.LogError("Failed to get OrbitEllipse from prefab:" + go.name);
            return;
        }
        // Awkward test code
        if (parent == star) {
            eb.paramBy = EllipseBase.ParamBy.AXIS_A;
            eb.a = Random.Range(minRadius, maxRadius);
            eb.inclination = Random.Range(-80f, 80f);
            eb.ecc = Random.Range(minEccentricity, maxEccentricity);
        } else {
            // moon test, so keep "a" small
            eb.paramBy = EllipseBase.ParamBy.AXIS_A;
            eb.a = moonRadius;
            eb.inclination = Random.Range(-80f, 80f);
            eb.ecc = 0;
        }
        if (bodyType == KEPLER) {
            eb.evolveMode = OrbitEllipse.evolveType.KEPLERS_EQN;
        }
        eb.Init();
        OrbitPredictor op = go.GetComponentInChildren<OrbitPredictor>();
        if (op != null) {
            op.centerBody = eb.centerObject;
        }
    }

    private void HyperInit(GameObject go, GameObject parent) {
        OrbitHyper hyper = go.GetComponent<OrbitHyper>();
        if (hyper == null) {
            Debug.LogError("Failed to get OrbitHyper from prefab:" + go.name);
            return;
        }
        // Awkward test code
        if (parent == star) {
            hyper.perihelion = Random.Range(minRadius, maxRadius);
            // aribitrary - start at fixed distance from peri
            hyper.r_initial = 1.0f * hyper.perihelion;
            hyper.inclination = Random.Range(-80f, 80f);
            hyper.ecc = Random.Range(1.1f, 2f);
            hyper.centerObject = parent;
        }
        hyper.Init();
        OrbitPredictor op = go.GetComponentInChildren<OrbitPredictor>();
        if (op != null) {
            op.centerBody = hyper.centerObject;
        }
    }

    private void AddBodyToParent(string bodyType, GameObject parent) {

        GameObject go;
        // getting long - RF later...
		if (bodyType == MASSLESS) {
            go = Instantiate(orbitingPrefab) as GameObject;
			masslessObjects.Add(go);
		} else if (bodyType == KEPLER) {
            go = Instantiate(orbitingPrefab) as GameObject;
			keplerObjects.Add(go);
        } else if (bodyType == BINARY) {
            go = Instantiate(binaryPrefab) as GameObject;
            binaryObjects.Add(go);
        } else if (bodyType == MASSIVE) {
            go = Instantiate(orbitingPrefab) as GameObject;
			massiveObjects.Add(go);
        } else if (bodyType == HYPER) {
            go = Instantiate(hyperPrefab) as GameObject;
            hyperObjects.Add(go);
        } else if (bodyType == MOON) {
            go = Instantiate(orbitingPrefab) as GameObject;
            moonObjects.Add(go);
        } else if (bodyType == FIXED) {
            go = Instantiate(fixedObject) as GameObject;
            fixedObjects.Add(go);
        } else {
            Debug.LogWarning("Do not understand string=" + bodyType);
            return;
        }

        NBody nbody = go.GetComponent<NBody>();
        if ((bodyType == MASSLESS) || (bodyType == MOON)) {
            nbody.mass = 0; 
        } else if (bodyType != FIXED) {
            nbody.mass = mass;
        }
        if (bodyType == FIXED) {
            // Add fixed stars in the X/Y plane at a fixed distance
            float theta = Random.Range(0, 2f*Mathf.PI);
            nbody.initialPos = new Vector3(fixedStarRadius * Mathf.Cos(theta),
                                    fixedStarRadius * Mathf.Sin(theta), 
                                    0);
        }
        go.transform.parent = parent.transform;

        if (bodyType == HYPER) {
            HyperInit(go, parent);

        } else if (bodyType != FIXED) {
            EllipseInit(go, parent, bodyType);
        }

        GravityEngine.instance.AddBody(go);

        TrailRenderer[] trails = go.GetComponentsInChildren<TrailRenderer>();
        foreach (TrailRenderer trail in trails) {
            trail.material.color = colors[colorIndex];
            colorIndex = (colorIndex + 1) % colors.Length;
            trail.enabled = true;
            trail.Clear();
        }

    }

     public void RemoveBody(string bodyType) {

		List<GameObject> bodyList = null;

        if (bodyType == MASSLESS) {
            bodyList = masslessObjects;
        } else if (bodyType == KEPLER) {
            bodyList = keplerObjects;
        } else if (bodyType == BINARY) {
            bodyList = binaryObjects;
        } else if (bodyType == HYPER) {
            bodyList = hyperObjects;
        } else if (bodyType == MOON) {
            bodyList = moonObjects;
        } else if (bodyType == FIXED) {
            bodyList = fixedObjects;
        } else {
			bodyList = massiveObjects;
		}
		if (bodyList.Count > 0) {
			int entry = (int)(Random.Range(0, (float) bodyList.Count));
			GameObject toDestroy = bodyList[entry];
            // massive body may have a moon, binary has two NBody kids
            if ((bodyType == BINARY) || (bodyType == MASSIVE)) {
                NBody[] bodies = toDestroy.GetComponentsInChildren<NBody>();
                foreach (NBody nb in bodies) {
                    GravityEngine.instance.RemoveBody(nb.gameObject);
                    if (bodyType == MASSIVE) {
                        moonObjects.Remove(nb.gameObject);
                    }
                }
            }
            if (bodyType != BINARY) { 
                GravityEngine.instance.RemoveBody(toDestroy);
            }
            bodyList.RemoveAt(entry);
			Destroy(toDestroy);
		} else {
			Debug.Log("All objects of that type removed.");
		}

	}

    void Update() {
        if (Input.GetKeyDown(KeyCode.C)) {
            GravityEngine.Instance().Clear();
            Debug.Log("Clear all bodies");
            List<GameObject>[] lists = { moonObjects,
                                        hyperObjects,
                                        fixedObjects,
                                        massiveObjects,
                                        masslessObjects,
                                        keplerObjects,
                                        binaryObjects};
            foreach (List<GameObject> l in lists) {
                foreach( GameObject g in l) {
                    Destroy(g);
                }
                l.Clear();
            }
        }
        if (runChaosMonkey) {
            RunChaosMonkey();
        }
    }

    /// <summary>
    /// Chaos Monkey to do random add/deletes at random time offsets. Goal it to make sure no errors
    /// or warning occur in the console. 
    /// </summary>

    private float nextAddTime = 0f;
    private void RunChaosMonkey() {
        if (Time.time > nextAddTime) {
            bool add = false;
            if (Random.Range(0f, 1f) < 0.5f) {
                add = true;
            }
            string type = adTypes[Random.Range(0, adTypes.Length)];
            if (add) {
                AddBody(type);
            } else {
                RemoveBody(type);
            }
            nextAddTime = Time.time + Random.Range(0f, 0.2f);
        }
    }

}
