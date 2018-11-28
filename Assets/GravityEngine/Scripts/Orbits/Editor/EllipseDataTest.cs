using UnityEngine;
using UnityEditor;
using NUnit.Framework;

public class EllipseDataTest {

	private static bool FloatEqual(float a, float b) {
		return (Mathf.Abs(a-b) < 1E-3); 
	}

	private static bool FloatEqual(float a, float b, double error) {
		return (Mathf.Abs(a-b) < error); 
	}

    [Test]
    // Create an NBody and check it's mass
    public void NbodyCreate()
    {
        const float mass = 1000f; 
        var gameObject = TestSetupUtils.CreateNBody(mass, new Vector3(0,0,0));

        NBody nbody = gameObject.GetComponent<NBody>();
        Assert.AreEqual(nbody.mass, nbody.mass);
    }

	[Test]
    // Does LH co-ordinate system do anything goofy to A x B ?
    // Seems not. 
    public void CrossProduct()
    {
    	Vector3 a = new Vector3(1f, 2f, 3f); 
    	Vector3 b = new Vector3(1f, 1f, 1f);
    	Vector3 c = Vector3.Cross(a,b);
        Assert.AreEqual(c, new Vector3(-1f, 2f, -1f));
    }

	[Test]
    // Create an NBody and check it's mass
    public void CircleA()
    {
        const float mass = 1000f; 
        GameObject star = TestSetupUtils.CreateNBody(mass, new Vector3(0,0,0));
        const float orbitRadius = 10f; 
		GameObject planet = TestSetupUtils.CreatePlanetInOrbit(star, 1f, orbitRadius);
        TestSetupUtils.SetupGravityEngine(star, planet);
        // confirm planet is in correct location
        Assert.AreEqual( Vector3.Distance(planet.transform.position, new Vector3(10f,0f,0)), 0);
		// take the velocity and check 
		OrbitData orbitData = new OrbitData();
		orbitData.SetOrbitForVelocity(planet.GetComponent<NBody>(), star.GetComponent<NBody>());
		Assert.AreEqual( orbitData.a, orbitRadius);
		Assert.AreEqual( orbitData.omega_uc, 0f);
    }

	[Test]
    // Check eccentricity and inclination
    public void CheckTestRV()
    {
        const float mass = 1000f; 
        GameObject star = TestSetupUtils.CreateNBody(mass, new Vector3(0,0,0));
        const float orbitRadius = 10f; 
		GameObject planet = TestSetupUtils.CreatePlanetInOrbit(star, 1f, orbitRadius);
        OrbitEllipse orbitEllipse = planet.GetComponent<OrbitEllipse>();
		orbitEllipse.ecc = .25f; 
		orbitEllipse.inclination = 25f;
		orbitEllipse.omega_uc = 10f;
		orbitEllipse.omega_lc = 20f;
		orbitEllipse.phase = 190f;

        OrbitData od = new OrbitData();
		od.a = orbitRadius; 
		od.ecc = 0.25f;
		od.inclination = 25f;
		od.omega_uc = 10f;
		od.omega_lc = 20f;
		od.phase = 190f;
        GravityEngine.Instance().UnitTestAwake();
		TestRV(od, planet, star, orbitRadius);
	}

	[Test]
    // Check eccentricity and inclination
    public void EllipseInclination()
    {
        const float mass = 1000f; 
        GameObject star = TestSetupUtils.CreateNBody(mass, new Vector3(0,0,0));
        const float orbitRadius = 10f; 
		GameObject planet = TestSetupUtils.CreatePlanetInOrbit(star, 1f, orbitRadius);
        OrbitEllipse orbitEllipse = planet.GetComponent<OrbitEllipse>();

		float eccentricity = 0.3f;
		orbitEllipse.ecc = eccentricity;
        TestSetupUtils.SetupGravityEngine(star, planet);

        // take the velocity and check 
        OrbitData orbitData = new OrbitData();
		orbitData.SetOrbitForVelocity(planet.GetComponent<NBody>(), star.GetComponent<NBody>());
		Assert.IsTrue( FloatEqual(orbitRadius, orbitData.a) );
		Assert.IsTrue( FloatEqual(eccentricity, orbitData.ecc) );

		// Try some values of inclination
		float[] inclinationValues = { 0f, 30f, 45f, 60f, 90f, 135f, 180f};
		foreach (float inc in inclinationValues) {
            Debug.Log("====EclipseInclination====    inc=" + inc);
			orbitEllipse.inclination = inc ;
            TestSetupUtils.SetupGravityEngine(star, planet);
            OrbitData od = new OrbitData();
			od.SetOrbitForVelocity(planet.GetComponent<NBody>(), star.GetComponent<NBody>());
			Debug.Log("TEST: incl = " + orbitEllipse.inclination + " ecc=" + orbitEllipse.ecc + " od:" + od.LogString());
            Assert.IsTrue( FloatEqual(inc, od.inclination) );
			TestRV(od, planet, star, orbitRadius);
		}
    }

	[Test]
    public void CirclePhaseNoInclination()
    {
        const float mass = 1000f; 
        GameObject star = TestSetupUtils.CreateNBody(mass, new Vector3(0,0,0));
        const float orbitRadius = 10f; 
		GameObject planet = TestSetupUtils.CreatePlanetInOrbit(star, 1f, orbitRadius);
		OrbitEllipse orbitEllipse = planet.GetComponent<OrbitEllipse>();

		// Try some values of om
		float[] phaseValues = { 30f, 45f, 60f, 90f, 135f, 180f, 210f, 340f};
		foreach (float phase in phaseValues) {
			orbitEllipse.phase = phase ; 
            TestSetupUtils.SetupGravityEngine(star, planet);
            OrbitData od = new OrbitData();
			od.SetOrbitForVelocity(planet.GetComponent<NBody>(), star.GetComponent<NBody>());
            orbitEllipse.Log("CirclePhaseNoInclination: ");
            Debug.Log("phase = " + phase + " od.phase=" + od.phase);
            Debug.Log(od.LogString());
            // Need a bit of leeway at 0 with error
            Assert.IsTrue( FloatEqual(phase, od.phase, 0.02) );
			TestRV(od, planet, star, orbitRadius);
		}
    }

	[Test]
    public void PhaseNoInclinationEllipse()
    {
        const float mass = 1000f; 
        GameObject star = TestSetupUtils.CreateNBody(mass, new Vector3(0,0,0));
        const float orbitRadius = 20f; 
		GameObject planet = TestSetupUtils.CreatePlanetInOrbit(star, 1f, orbitRadius);
		OrbitEllipse orbitEllipse = planet.GetComponent<OrbitEllipse>();
		orbitEllipse.ecc = 0.2f;

        // Try some values of om
        float[] phaseValues = { 30f, 45f, 60f, 90f, 135f, 180f, 225f, 270f, 325f,  0f};
		foreach (float phase in phaseValues) {
			orbitEllipse.phase = phase ;
            TestSetupUtils.SetupGravityEngine(star, planet);

            OrbitData od = new OrbitData();
			od.SetOrbitForVelocity(planet.GetComponent<NBody>(), star.GetComponent<NBody>());
			Debug.Log("\n***Test Case: phase = " + phase + "\n OD:=" + od.LogString() );
			// Need a bit of leeway at 0 with error
			Assert.IsTrue( FloatEqual(phase, od.phase, 0.02) );
            Assert.IsTrue( FloatEqual(0f, od.omega_lc, 0.02));
            Assert.IsTrue( FloatEqual(0f, od.omega_uc, 0.02));
            TestRV(od, planet, star, orbitRadius);
		}
    }

	[Test]
    public void PhaseInclinedEllipse()
    {
        const float mass = 1000f; 
        GameObject star = TestSetupUtils.CreateNBody(mass, new Vector3(0,0,0));
        const float orbitRadius = 10f; 
		GameObject planet = TestSetupUtils.CreatePlanetInOrbit(star, 1f, orbitRadius);
        OrbitEllipse orbitEllipse = planet.GetComponent<OrbitEllipse>();
		orbitEllipse.ecc = 0.4f;
		orbitEllipse.inclination = 60f;

		// Try some values of om
		float[] phaseValues = { 30f, 45f, 60f, 90f, 135f, 180f, 0f};
		foreach (float phase in phaseValues) {
			orbitEllipse.phase = phase ;
            TestSetupUtils.SetupGravityEngine(star, planet);
            OrbitData od = new OrbitData();
			od.SetOrbitForVelocity(planet.GetComponent<NBody>(), star.GetComponent<NBody>());
			Debug.Log("phase = " + phase + " od.phase=" + od.phase);
			// Need a bit of leeway at 0 with error
			Assert.IsTrue( FloatEqual(phase, od.phase, 0.02) );
			TestRV(od, planet, star, orbitRadius);
		}
    }

	[Test]
    public void PhaseRetrogradeEllipse()
    {
        const float mass = 1000f; 
        GameObject star = TestSetupUtils.CreateNBody(mass, new Vector3(0,0,0));
        const float orbitRadius = 10f; 
		GameObject planet = TestSetupUtils.CreatePlanetInOrbit(star, 1f, orbitRadius);
        OrbitEllipse orbitEllipse = planet.GetComponent<OrbitEllipse>();
		orbitEllipse.ecc = 0.4f;
		orbitEllipse.inclination = 180f;

		// Try some values of om
		float[] phaseValues = { 30f, 45f, 60f, 90f, 135f, 180f, 0f};
		foreach (float phase in phaseValues) {
			orbitEllipse.phase = phase ;
            TestSetupUtils.SetupGravityEngine(star, planet);
            OrbitData od = new OrbitData();
			od.SetOrbitForVelocity(planet.GetComponent<NBody>(), star.GetComponent<NBody>());
			Debug.Log("phase = " + phase + " od.phase=" + od.phase);
			// Need a bit of leeway at 0 with error
			Assert.IsTrue( FloatEqual(phase, od.phase, 0.02) );
			TestRV(od, planet, star, orbitRadius);
		}
    }

	[Test]
    public void CirclePhaseOmegaInclined()
    {
        const float mass = 1000f; 
        GameObject star = TestSetupUtils.CreateNBody(mass, new Vector3(0,0,0));
        const float orbitRadius = 10f; 
		GameObject planet = TestSetupUtils.CreatePlanetInOrbit(star, 1f, orbitRadius);
        OrbitEllipse orbitEllipse = planet.GetComponent<OrbitEllipse>();
		orbitEllipse.inclination = 20f;

		// Try some values of phase
		float[] phaseValues = { 0f, 30f, 45f, 60f, 90f, 135f, 180f, 210f, 320f};
		float[] omegaUValues = { 0f, 30f, 45f, 60f, 90f, 135f, 180f, 210f, 320f};
		foreach (float phase in phaseValues) {
			foreach (float omegau in omegaUValues) {
				orbitEllipse.phase = phase ; 
				orbitEllipse.omega_uc = omegau;
                TestSetupUtils.SetupGravityEngine(star, planet);
                orbitEllipse.Log("Initial circle:");
				OrbitData od = new OrbitData();
				od.SetOrbitForVelocity(planet.GetComponent<NBody>(), star.GetComponent<NBody>());
				Debug.Log("TEST: phase = " + phase + " Omega=" + omegau + " od:" + od.LogString() );
				TestRV(od, planet, star, orbitRadius);
			}
		}
    }

	// Do omega_uc with NO inclination. This will come back as omega_lc, since
    // Omega is angle to line of nodes - and does not apply if inclination=0
	[Test]
    public void OmegaUNoInclination()
    {
        const float mass = 1000f; 
        GameObject star = TestSetupUtils.CreateNBody(mass, new Vector3(0,0,0));
        const float orbitRadius = 10f; 
		GameObject planet = TestSetupUtils.CreatePlanetInOrbit(star, 1f, orbitRadius);
        OrbitEllipse orbitEllipse = planet.GetComponent<OrbitEllipse>();
		orbitEllipse.ecc = 0.1f;

		// Try some values of om
		float[] omegaValues = { 0f, 30f, 45f, 60f, 90f, 135f, 180f, 210f, 320f};
		foreach (float omega in omegaValues) {
			orbitEllipse.omega_uc = omega ;
            TestSetupUtils.SetupGravityEngine(star, planet);

            OrbitData od = new OrbitData();
			od.SetOrbitForVelocity(planet.GetComponent<NBody>(), star.GetComponent<NBody>());
			Debug.Log("Omega = " + omega + " od.omega_lc=" + od.omega_lc + " od:" + od.LogString());
			Assert.IsTrue( FloatEqual(omega, od.omega_lc, 0.4) );
		}
    }
	// Do omega_uc with inclination
	// Result may be a different (but equivilent) choice of Euler angles for omegas
	[Test]
    public void OmegaUEllipseInclined()
    {
        const float mass = 1000f; 
        GameObject star = TestSetupUtils.CreateNBody(mass, new Vector3(0,0,0));
        const float orbitRadius = 10f; 
		GameObject planet = TestSetupUtils.CreatePlanetInOrbit(star, 1f, orbitRadius);
        OrbitEllipse orbitEllipse = planet.GetComponent<OrbitEllipse>();
		orbitEllipse.InitNBody(1f, 1f);

		orbitEllipse.inclination = 30f;
		orbitEllipse.ecc = 0.2f;
		orbitEllipse.Init();
		orbitEllipse.ApplyScale(1f);
		orbitEllipse.InitNBody(1f, 1f);

		// Try some values of om
		float[] omegaValues = { 0f, 30f, 45f, 60f, 90f, 135f, 180f, 210f, 310f};
		foreach (float omega in omegaValues) {
			orbitEllipse.omega_uc = omega ;
            TestSetupUtils.SetupGravityEngine(star, planet);

            orbitEllipse.Log("OmegaU inclined:");
			OrbitData od = new OrbitData();
			od.SetOrbitForVelocity(planet.GetComponent<NBody>(), star.GetComponent<NBody>());
			Debug.Log("omegaU = " + omega + " od.omega_uc=" + od.omega_uc);
			TestRV(od, planet, star, orbitRadius);
		}
    }

	// Do omega_lc with inclination
	// Result may be a different (but equivilent) choice of Euler angles for omegas
	[Test]
    public void OmegaLEllipseInclined()
    {
        const float mass = 1000f; 
        GameObject star = TestSetupUtils.CreateNBody(mass, new Vector3(0,0,0));
        const float orbitRadius = 10f; 
		GameObject planet = TestSetupUtils.CreatePlanetInOrbit(star, 1f, orbitRadius);
        OrbitEllipse orbitEllipse = planet.GetComponent<OrbitEllipse>();
		orbitEllipse.InitNBody(1f, 1f);

		orbitEllipse.inclination = 30f;
		orbitEllipse.ecc = 0.2f;
		orbitEllipse.ApplyScale(1f);
		orbitEllipse.Init();
		orbitEllipse.InitNBody(1f, 1f);

		// Try some values of om
		float[] omegaValues = { 0f, 30f, 45f, 60f, 90f, 135f, 180f, 210f, 310f};
		foreach (float omega in omegaValues) {
			orbitEllipse.omega_lc = omega ;
            TestSetupUtils.SetupGravityEngine(star, planet);

            orbitEllipse.Log("\n*** Case: OmegaL inclined: " + omega + "\n");
			OrbitData od = new OrbitData();
			od.SetOrbitForVelocity(planet.GetComponent<NBody>(), star.GetComponent<NBody>());
			Debug.Log("omegaL = " + omega + " od.omega_lc=" + od.omega_lc + "\n OD:" + od.LogString() );
			TestRV(od, planet, star, orbitRadius);
		}
    }

	[Test]
    public void OmegasEllipseInclination()
    {
        const float mass = 1000f; 
        GameObject star = TestSetupUtils.CreateNBody(mass, new Vector3(0,0,0));
        const float orbitRadius = 10f; 
		GameObject planet = TestSetupUtils.CreatePlanetInOrbit(star, 1f, orbitRadius);
        OrbitEllipse orbitEllipse = planet.GetComponent<OrbitEllipse>();
		orbitEllipse.InitNBody(1f, 1f);

		orbitEllipse.inclination = 35f;
		orbitEllipse.ecc = 0.25f;
		orbitEllipse.Init();
		orbitEllipse.ApplyScale(1f);
		orbitEllipse.InitNBody(1f, 1f);

		// Try some values of om
		float[] omegaValues = { 0f, 30f, 45f, 60f, 90f, 135f, 180f};
		foreach (float omega in omegaValues) {
			orbitEllipse.omega_uc = omega ;
            TestSetupUtils.SetupGravityEngine(star, planet);

            OrbitData od = new OrbitData();
			od.SetOrbitForVelocity(planet.GetComponent<NBody>(), star.GetComponent<NBody>());
			TestRV(od, planet, star, orbitRadius);
		}
    }

    [Test]
    public void OmegaULEllipseInclined() {
        const float mass = 1000f;
        GameObject star = TestSetupUtils.CreateNBody(mass, new Vector3(0, 0, 0));
        const float orbitRadius = 10f;
        GameObject planet = TestSetupUtils.CreatePlanetInOrbit(star, 1f, orbitRadius);
        OrbitEllipse orbitEllipse = planet.GetComponent<OrbitEllipse>();
        orbitEllipse.InitNBody(1f, 1f);

        orbitEllipse.inclination = 30f;
        orbitEllipse.ecc = 0.2f;
        orbitEllipse.Init();
        orbitEllipse.ApplyScale(1f);
        orbitEllipse.InitNBody(1f, 1f);

        // Try some values of om
        float[] omegaUValues = { 0f, 30f, 45f, 60f, 90f, 135f, 180f, 210f, 310f };
        float[] omegaLValues = { 0f, 30f, 45f, 60f, 90f, 135f, 180f, 210f, 310f };
        foreach (float omegaU in omegaUValues) {
            foreach (float omegaL in omegaLValues) {
                orbitEllipse.omega_uc = omegaU;
                orbitEllipse.omega_lc = omegaL;
                TestSetupUtils.SetupGravityEngine(star, planet);

                orbitEllipse.Log("OmegaUL inclined:");
                OrbitData od = new OrbitData();
                od.SetOrbitForVelocity(planet.GetComponent<NBody>(), star.GetComponent<NBody>());
                Debug.Log("omegaU = " + omegaU + "omegaL = " + omegaL + " OD:=" + od.LogString() );
                TestRV(od, planet, star, orbitRadius);
            }
        }
    }

    private void TestRV(OrbitData od, GameObject planet, GameObject star, float orbitRadius)  {
 
		GameObject testPlanet = TestSetupUtils.CreatePlanetInOrbit(star, 1f, orbitRadius);
        OrbitEllipse testEllipse = testPlanet.GetComponent<OrbitEllipse>();
        // Run init explicitly to update transform details
		testEllipse.InitFromOrbitData( od);
 
        // Awkward but cannot add a new object to GE when it is stopped, so re-add all three
        GravityEngine ge = GravityEngine.Instance();
        ge.Clear();
        ge.AddBody(star);
        ge.AddBody(planet);
        ge.AddBody(testPlanet);
        ge.Setup();
        ge.LogDump();
        Vector3 r_od = ge.GetPhysicsPosition(testPlanet.GetComponent<NBody>());
		Vector3 v_od = ge.GetVelocity(testPlanet);
        Vector3 r_i = ge.GetPhysicsPosition(planet.GetComponent<NBody>());
        Vector3 v_i = ge.GetVelocity(planet);
        Debug.Log(" r_i=" + r_i + " r_od=" + r_od + " delta=" + Vector3.Distance(r_i, r_od));
		Debug.Log(" v_i=" + v_i + " v_od=" + v_od + " delta=" + Vector3.Distance(v_i, v_od));
		Assert.IsTrue( FloatEqual( Vector3.Distance(r_i, r_od), 0f, 1E-2));
		Assert.IsTrue( FloatEqual( Vector3.Distance(v_i, v_od), 0f, 1E-2));

    }

}
