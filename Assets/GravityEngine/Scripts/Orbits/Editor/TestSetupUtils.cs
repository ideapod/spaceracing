using UnityEngine;
using System.Collections;

public class TestSetupUtils : MonoBehaviour {

	public static GameObject CreateNBody(float mass, Vector3 position) {
		GameObject gameObject = new GameObject();
		gameObject.transform.position = position;

		NBody nbody = gameObject.AddComponent<NBody>();
		nbody.mass = mass;
		return gameObject;
	}

	// Create a planet in orbit around center object with semi-major axis a
	public static GameObject CreatePlanetInOrbit(GameObject center, float mass, float a) {
		// position will be trumped by orbit
		GameObject planet = CreateNBody(mass, new Vector3(1,0,0));

		OrbitEllipse orbitEllipse = planet.AddComponent<OrbitEllipse>();
		orbitEllipse.a = a;
		orbitEllipse.SetCenterBody(center);
        return planet;
	}

	// Create a planet in orbit around center object with semi-major axis a
	public static GameObject CreatePlanetInHyper(GameObject center, float mass) {
		// position will be trumped by orbit
		GameObject planet = CreateNBody(mass, new Vector3(1,0,0));

		OrbitHyper orbitHyper = planet.AddComponent<OrbitHyper>();
		orbitHyper.SetCenterBody(center);
		return planet;
	}

    public static void SetupGravityEngine(GameObject centerBody, GameObject orbitingBody) {
        GravityEngine ge = GravityEngine.Instance();
        if (ge == null)
            Debug.LogError("No GE in scene");
        if (ge.evolveAtStart) {
            Debug.LogError("Evolve at start set. Are you in the TestRunner scene?");
        } else if (ge.detectNbodies) {
            Debug.LogError("Detect NBodies at start set. Are you in the TestRunner scene?");
        }
        ge.UnitTestAwake();
        ge.Clear();
        ge.AddBody(centerBody);
        ge.AddBody(orbitingBody);
        ge.Setup();
        ge.LogDump();
    }
}
