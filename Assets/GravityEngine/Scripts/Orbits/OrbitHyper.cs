using UnityEngine;
#if UNITY_EDITOR
using UnityEditor;
#endif
using System.Collections;

/// <summary>
///
/// Class used to handle the hyperbolic orbital parameters and draw the orbit path in the editor using Gizmos. 
///
/// How to specify an hyperbola in 3D:
///
///    p - pericenter - distance from focus of hyperbola to point of closest approach to that focus
///
/// shape: controlled by ecc (eccentricity) 0 for a circle, 0.99 for a very long thin ellipse
///
/// orientation:
///  The standard orbit parameters are used. You can develop some intuition for these by chaging them in the Editor and observing
///  the change in the orbit. 
///
///  Orientation is defined with respect to the positive X axis. 
///    omega (lower case) - is a rotation in the plane of the orbit 
///    Inclination - is the tilt of the closest approach vector to the XY plance
///    Omega (capital Omega) - is the rotation around Z after preceeding rotations
///
/// </summary>
public class OrbitHyper : MonoBehaviour, INbodyInit, IOrbitPositions, IOrbitScalable  {

	public enum evolveType {GRAVITY_ENGINE, KEPLERS_EQN};
	//! Use GRAVITY_ENGINE to evolve or move in a fixed KEPLER orbit. 
	public evolveType evolveMode = evolveType.GRAVITY_ENGINE;

	//! object to orbit around (if null, will take parent game object)
	public GameObject centerObject; 
	
	// Orbit parameters (user control via FixedEllipseEditor)
	// These parameters are in world space. 
	//! eccentricity (0..1, 0=circle, 1=linear)
	public float ecc = 2.0f; 			

	/// <summary>
	/// Hyperbola parameters:
	/// The definition is typically in terms of peri[helion] (closest approach to e.g. Sun) but the
	/// hyperbola equation uses a.
	///
	/// a is calculated from the perihelion and used for orbital calculations.
	///
	/// </summary>

	//! point of closest approach
	public float perihelion = 10f; 		
	//! point of closest approach
	private float perihelion_scaled = 10f;

    //! fraction of the branch of the hyperbola to display in OrbitPositions
    public float branchDisplayFactor = 0.5f;

	//! semi-major axis - based on paramBy user can specify a OR p. a = p/(1-ecc)
	private float a_scaled; 			

	//! "longitude of ascending node" - angle from x-axis to line from focus to pericenter
	public float omega_uc; 		
	//! "argument of perienter" - angle from ascending node to pericenter
	public float omega_lc; 		
	//! inclination (degrees!)
	public float inclination; 	
	//! initial distance from focus
	public float r_initial = 10f;
    private float r_initial_scaled = 10f; 

	//! initial distance on outbound leg
	public bool r_initial_outbound = false; 	

	private float b; // hyperbola semi-minor axis

    private Vector3 centerPos;

	protected Quaternion hyper_orientation;

	protected Vector3 xunit = Vector3.right;
	protected Vector3 yunit = Vector3.up;
	protected Vector3 zunit = Vector3.forward;

	protected NBody centerNbody;

    void Start() {
		Init();
	}
	
	/// <summary>
	/// Init the hyperbola, verify a center body is present and determine orientation.
	/// </summary>
	public void Init () {
        CalcOrbitParams();
        centerNbody = OrbitUtils.GetCenterNbody(transform, centerObject);
		CalculateRotation();

        NBody nbody = GetComponent<NBody>();
        // particle ring would not have an NBody
        if (nbody != null) {
            SetInitialPosition(nbody);
        }
    }

    public void InitFromOrbitData(OrbitData od) {
		a_scaled = od.a; 
		perihelion_scaled = od.perihelion;  
		ecc = od.ecc; 
		omega_lc = od.omega_lc;
		omega_uc = od.omega_uc; 
		inclination = od.inclination;
		r_initial = od.r_initial;
		Init();
	}

	/// <summary>
	/// Sets the center body and initializes the hyperbola configuration.
	/// </summary>
	/// <param name="centerBody">Center body.</param>
	public void SetCenterBody(GameObject centerBody) {
		centerObject = centerBody; 
		Init(); 
	}
	
	// Update the derived orbit parameters
	protected void CalcOrbitParams() {
		// Protect against scripts changing ecc to < 1
		if (ecc < 1f) {
		    Debug.LogWarning("Detected ecc < 1. Set to 1.01");
		    ecc = 1.01f; 
		}
		// Roy has this as p = a(e^2-1) but says when f=0, r=a(e-1)
		// Wolfram has p = a (e^2-1)/e
        // p is semi-latus rectum, not periapse
		// p is the distance at f = +/- Pi/2  i.e. where the line x=const through the focus
		// intercepts the hyperbola. The point of closest approach is when f=0, cosf=1, r=perihelion/(1+e)
		a_scaled = perihelion_scaled/(ecc-1f);
		b = Mathf.Sqrt( a_scaled*a_scaled * (ecc*ecc-1f));
	}

	protected void CalculateRotation() {
		// Following Murray and Dermot Ch 2.8 Fig 2.14
		// Quaternions go L to R (matrices are R to L)
		hyper_orientation = Quaternion.AngleAxis(omega_uc, zunit ) *
							  Quaternion.AngleAxis(inclination, xunit) * 
							  Quaternion.AngleAxis(omega_lc, zunit);
	}



    private void SetInitialPosition(NBody nbody) {
		// Phase is TRUE anomoly f
		float f = ThetaForR(r_initial_scaled);
		// Murray and Dermot (2.20)
		float r = r_initial_scaled;
		Vector3 position_xy = new Vector3( r * Mathf.Cos(f), r* Mathf.Sin(f), 0);
        // move from XY plane to the orbital plane and scale to world space
        position_xy = hyper_orientation * position_xy;

        // orbit position is WRT center
        // orbit position is WRT center. Could be adding dynamically to an object in motion, so need current position. 
        Vector3 centerPos = Vector3.zero;
        // used by widgets - so need to get explcitly
        centerNbody = OrbitUtils.GetCenterNbody(this.transform, centerObject);
        if (centerNbody.engineRef != null) {
            centerPos = GravityEngine.Instance().GetPhysicsPosition(centerNbody);
        } else {
            // setup - not yet added to GE
            centerPos = centerNbody.initialPhysPosition;
        }
        nbody.initialPhysPosition = position_xy + centerPos;

    }

    private float ThetaForR(float r) {
		float theta = 0; 
		// solve r_i = a(e^2-1)/(1+eCos(theta)) for theta
		if (Mathf.Abs(r_initial_scaled) > 1E-2) {
			float arg = ((a_scaled*(ecc*ecc-1))/ r_initial_scaled - 1f)/ecc;
			arg = Mathf.Max(-1f, arg);
			arg = Mathf.Min(1f, arg);
			theta = Mathf.Acos(arg);
		}
		return theta;
	}

	private float RforTheta(float theta) {
		return a_scaled * ( ecc* ecc - 1)/(1f + ecc * Mathf.Cos(theta));
	}

	public bool IsFixed() {
		return (evolveMode == evolveType.KEPLERS_EQN);
	}

	/// <summary>
	/// Inits the N body position and velocity based on the hyperbola parameters and the 
	/// position and velocity of the parent. 
	/// </summary>
	/// <param name="physicalScale">Physical scale.</param>
	public void InitNBody(float physicalScale, float massScale) {

        Init();
		float a_phy = a_scaled/physicalScale;
		NBody nbody = GetComponent<NBody>();
		
		// Phase is TRUE anomoly f
		float f = ThetaForR(r_initial_scaled);
		// Murray and Dermot (2.20)
		float n = Mathf.Sqrt( (float)(centerNbody.mass * massScale)/(a_phy*a_phy*a_phy));
		float denom = Mathf.Sqrt( ecc*ecc - 1f);
		// reverse sign from text to get prograde motion
		float xdot = 1f * n * a_phy * Mathf.Sin(f)/denom;
		float ydot = -1f * n * a_phy * (ecc + Mathf.Cos(f))/denom;
		if (!r_initial_outbound) {
			xdot *= -1f;
			ydot *= -1f;
		}

		// Init functions are called in the engine by SetupOneBody and calls of parent vs children/grandchildren etc.
		// can be in arbitrary order. A child will need info from parent for position and velocity. Ensure parent
		// has inited.
		INbodyInit centerInit = centerObject.GetComponent<INbodyInit>();
		if (centerInit != null) {
			centerInit.InitNBody(physicalScale, massScale);
		}
        if (centerNbody.engineRef != null) {
            centerPos = GravityEngine.Instance().GetPhysicsPosition(centerNbody);
        } else {
            // setup - not yet added to GE
            centerPos = centerNbody.initialPhysPosition;
        }

        Vector3 v_xy = new Vector3( xdot, ydot, 0);
		Vector3 vphy = hyper_orientation * v_xy + centerNbody.vel_phys;
		nbody.vel_phys = vphy;
	}	

	private Vector3 PositionForTheta(float theta) {
		float r = RforTheta(theta);
		Vector3 position = new Vector3( -r * Mathf.Cos (theta), r * Mathf.Sin (theta), 0);
		// move from XY plane to the orbital plane
		Vector3 newPosition = hyper_orientation * position; 
		// orbit position is WRT center
		newPosition += centerPos;
        return newPosition;
	}

	private Vector3 PositionForThetaLeftBranch(float theta, Vector3 centerPos) {
		float r = RforTheta(theta);
		// flip sign of X to get left branch
		Vector3 position = new Vector3( r * Mathf.Cos (theta), r * Mathf.Sin (theta), 0);
		// move from XY plane to the orbital plane
		Vector3 newPosition = hyper_orientation * position; 
		// orbit position is WRT center
		newPosition += centerPos;
		return newPosition;
	}

	private Vector3 PositionForY(float y) {
		float x = a_scaled * Mathf.Sqrt( 1 + y*y/(b*b));
		// focus is at x = -(a*e), want to translate to origin is at focus
		Vector3 position = new Vector3( -x + a_scaled*ecc, y, 0);
		// move from XY plane to the orbital plane
		Vector3 newPosition = hyper_orientation * position; 
		// orbit position is WRT center
		newPosition += centerPos;
        return newPosition;
	}

	/// <summary>
	/// Calculate an array of orbit positions. Used by the OrbitPredictor, OrbitRenderer and Editor
	/// Gimzo to illustrate the hyperbola. 
	/// </summary>
	/// <returns>The positions.</returns>
	/// <param name="numPoints">Number points.</param>
	public Vector3[] OrbitPositions(int numPoints, Vector3 centerPos, bool doSceneMapping) {

        CalculateRotation();

        Vector3[] emptyArray = {new Vector3(0,0,0), new Vector3(0,0,0)};
		// need to have a center to create positions.
		if (centerObject == null) {
			centerObject = transform.parent.gameObject;
			if (centerObject == null) {
				return emptyArray;
			}
		}
		Vector3[] points = new Vector3[numPoints];
		float theta = -1f*branchDisplayFactor * Mathf.PI;
        float dTheta = 2f* Mathf.Abs(theta) / (float)numPoints;
        GravityEngine ge = GravityEngine.Instance();
		for (int i=0; i < numPoints; i++)
		{
			points[i] = PositionForThetaLeftBranch(theta, centerPos);
            if (NUtils.VectorNaN(points[i])) {
                points[i] = Vector3.zero;
            } else  if (doSceneMapping && ge.mapToScene) {
                points[i] = ge.MapToScene(points[i]);
            }
			theta += dTheta;
		}
		return points;
	}

	public void Log(string prefix) {
		Debug.Log(string.Format("orbitHyper: {0} a_scaled={1} ecc={2} peri={3} i={4} Omega={5} omega={6} r_initial={7}", 
								prefix, a_scaled, ecc, perihelion, inclination, omega_uc, omega_lc, r_initial));
	}

	/// <summary>
	/// Apply scale to the orbit. This is used by the inspector scripts during
	/// scene setup. Do not use at run-time.
	/// </summary>
	/// <param name="scale">Scale.</param>
	public void ApplyScale(float scale) {
		perihelion_scaled = perihelion * scale;
        r_initial_scaled = r_initial * scale;
		CalcOrbitParams();
        NBody nbody = GetComponent<NBody>();
        SetInitialPosition(nbody);
	}

    /// <summary>
    /// Return the center object around which this ellipse is defined.
    /// </summary>
    /// <returns>The center object.</returns>
    public GameObject GetCenterObject() {
        // need to have a center to draw gizmo.
        if (centerObject == null && transform.parent != null) {
            centerObject = transform.parent.gameObject;
        }
        centerNbody = centerObject.GetComponent<NBody>();
        if (centerNbody == null) {
            Debug.LogError("centerBody does not have NBody component: " + centerObject.name);
        }
        return centerObject;
    }

#if UNITY_EDITOR
    /// <summary>
    /// Displays the path of the elliptical orbit when the object is selected in the editor. 
    /// </summary>
    void OnDrawGizmosSelected()
	{
        // need to have a center to draw gizmo.
        if (GetCenterObject() == null) {
            return;
        }
        // need to have a center to draw gizmo.
        centerNbody = centerObject.GetComponent<NBody>();
        if (centerNbody == null) {
            return;
        }
        // only display if this object is directly selected
        if (Selection.activeGameObject != transform.gameObject) {
			return;
		}
		CalcOrbitParams();
		CalculateRotation();

		const int NUM_STEPS = 100; 
		const int STEPS_PER_RAY = 20; 
		int rayCount = 0; 
		Gizmos.color = Color.white;

        // Center object may need to determine it's position in an orbit
        // and update it's intialPhyPosition
        GravityEngine ge = GravityEngine.Instance();
        centerNbody.InitPosition(ge);
        centerNbody.EditorUpdate(ge);
        Vector3 centerPos = centerNbody.transform.position;
        Init();
        // use transform (can't ask GE it's not setup yet)
        Vector3[] points = OrbitPositions(NUM_STEPS, centerPos, false);

		for (int i=1; i < NUM_STEPS; i++) {
			Gizmos.DrawLine(points[i-1], points[i] );
			// draw rays from focus
			rayCount = (rayCount+1)%STEPS_PER_RAY;
			if (rayCount == 0) {
				Gizmos.DrawLine(centerPos, points[i] );
			}
		}
		Gizmos.color = Color.white;
		// Draw the axes in a different color
		Gizmos.color = Color.red;
		Gizmos.DrawLine(PositionForTheta(0.5f*Mathf.PI), PositionForTheta(-0.5f*Mathf.PI) );
		Gizmos.color = Color.blue;
		Gizmos.DrawLine(PositionForY(0f), centerObject.transform.position );

        // move body to location specified by parameters
        if (!ge.enabled) {
            NBody nbody = GetComponent<NBody>();
            if (nbody != null) {
                nbody.EditorUpdate(ge);
            }
        }

    }

#endif	
}

