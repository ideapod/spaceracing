using UnityEngine;
using System.Collections;

/// <summary>
/// Calculate the correct initial conditions for a specified elliptical orbit around a designated body OR
/// moves the object in a fixed ellipse determined by Kepler's equation. 
///
/// In Gravity Engine mode the orbit will evolve based on the global gravitation field and due to interactions with other bodies
/// the actual orbit may not be the ellipse shown in the Unity editor due to these perturbations. 
///
/// In KEPLER mode the orbit will be constrained to move on the ellipse specified. 
///
/// This script must be attached to a gameObject with an NBody component. 
///
/// </summary>

[RequireComponent(typeof(NBody))]
public class EllipseDouble : MonoBehaviour, INbodyInit, IFixedOrbit, IOrbitScalable  {

    public enum ParamBy { AXIS_A, CLOSEST_P };

    //! Define ellipse by semi0major axis (A) or closest approach (P)
    public ParamBy paramBy = ParamBy.AXIS_A;

    public enum evolveType {GRAVITY_ENGINE, KEPLERS_EQN};

	//! Use GRAVITY_ENGINE to evolve or move in a fixed KEPLER orbit. 
	public evolveType evolveMode = evolveType.GRAVITY_ENGINE;

    public double a;
    public double a_scaled;
    public double p;
    public double p_scaled;
    public double omega_uc;
    public double omega_lc;

    public double ecc;
    public double inclination;
    public double phase;

    public GameObject centerObject; 

    //! Current world position
    private Vector3 position; 

    //! Current velocity in scaled units
    private Vector3 velocityScaled;
    private NBody nbody;

    private NBody centerNbody; 

	public bool IsFixed() {
		return (evolveMode == evolveType.KEPLERS_EQN);
	}

    private Quaternion ellipse_orientation = Quaternion.identity;

    // ODD. Mathf has this but not Math
    private static double DEG_TO_RAD = System.Math.PI / 180.0;

    private void Init() {

        // mode determines if user wants to define wrt to A or P. Determine other quantity 
        if (paramBy == ParamBy.AXIS_A) {
            p = a * (1 - ecc);
        } else if (paramBy == ParamBy.CLOSEST_P) {
            a = p / (1 - ecc);
        }
        centerNbody = OrbitUtils.GetCenterNbody(transform, centerObject);
        CalculateRotation();
        NBody nbody = GetComponent<NBody>();
        // particle ring would not have an NBody
        if (nbody != null) {
            SetInitialPosition(nbody);
        }
    }

    private void CalculateRotation() {
        // TODO - double version of Quaternions
    }

    /// <summary>
    /// Inits the N body position and velocity based on the ellipse parameters and the 
    /// position and velocity of the parent. 
    /// </summary>
    /// <param name="physicalScale">Physical scale.</param>
    public void InitNBody(float physicalScale, float massScale) {

        centerNbody = OrbitUtils.GetCenterNbody(this.transform, centerObject);
		double a_phy = a/physicalScale;
		nbody = GetComponent<NBody>();

        // Phase is TRUE anomoly f
        double f = phase * DEG_TO_RAD; 
        // Murray and Dermot 
        // (2.26)
        // This should really be (M+m), but assume m << M
        // (massScale is added in at the GE level)
        double n = System.Math.Sqrt( (centerNbody.mass * massScale)/(a_phy*a_phy*a_phy));
		// (2.36)
		double denom = System.Math.Sqrt( 1.0 - ecc*ecc);
        double xdot = -1f * n * a_phy * System.Math.Sin(f)/denom;
        double ydot = n * a_phy * (ecc + System.Math.Cos(f))/denom;

		// Init functions are called in the engine by SetupOneBody and calls of parent vs children/grandchildren etc.
		// can be in arbitrary order. A child will need info from parent for position and velocity. Ensure parent
		// has inited. Init may get called more than once. 
		INbodyInit centerInit = centerObject.GetComponent<INbodyInit>();
		if (centerInit != null) {
			centerInit.InitNBody(physicalScale, massScale);
		}

		Vector3 v_xy = new Vector3( (float) xdot, (float) ydot, 0);
        // if we're evolving, get the latest velocity
        if (centerNbody.engineRef != null) {
            centerNbody.UpdateVelocity();
        }
		Vector3 vphy = ellipse_orientation * v_xy + centerNbody.vel_phys;
		nbody.vel_phys = vphy;
		SetInitialPosition(nbody);
	}

    /// <summary>
    /// Sets the initial position based on the orbit parameters. Used in the init phase to set the NBody in the
    /// correct position in the scene before handing control GE. 
    /// </summary>
    private void SetInitialPosition(NBody nbody) {

        double phaseRad = phase * DEG_TO_RAD;
        // position object using true anomoly (angle from  focus)
        double r = a_scaled * (1f - ecc * ecc) / (1f + ecc * System.Math.Cos(phaseRad));

        Vector3 pos = new Vector3((float)(r * System.Math.Cos(phaseRad)),
                                  (float)(r * System.Math.Sin(phaseRad)), 
                                  0);
        // move from XY plane to the orbital plane
        Vector3 new_p = ellipse_orientation * pos;
        // orbit position is WRT center. Could be adding dynamically to an object in motion, so need current position. 
        Vector3 centerPos = Vector3.zero;
        // used by widgets - so need to get explcitly
        centerNbody = OrbitUtils.GetCenterNbody(transform, centerObject);
        if (centerNbody.engineRef != null) {
            centerPos = GravityEngine.Instance().GetPhysicsPosition(centerNbody);
        } else {
            // setup - not yet added to GE
            centerPos = centerNbody.initialPhysPosition;
        }
        nbody.initialPhysPosition = new_p + centerPos;

    }

    /// <summary>
    /// Inits from solar body.
    /// </summary>
    /// <param name="sbody">Sbody.</param>
    //public void InitFromSolarBody(SolarBody sbody) {
    //	a = sbody.a; 
    //	ecc = sbody.ecc; 
    //	omega_lc = sbody.omega_lc;
    //	omega_uc = sbody.omega_uc; 
    //	inclination = sbody.inclination;
    //	phase = sbody.longitude;
    //	Init();
    //	ApplyScale(GravityEngine.Instance().GetLengthScale());
    //}

    // Fixed motion code

    private double a_phy; 	// semi-major axis scaled for physics space
	private double orbitPeriod; 
	private double mean_anomoly_phase;
	private GravityEngine gravityEngine;

    public double GetPeriod() {
        if (centerNbody == null)
            return 0;
        // Use Find to allow Editor to use method. 
        if (gravityEngine == null) {
            gravityEngine = (GravityEngine)FindObjectOfType(typeof(GravityEngine));
            if (gravityEngine == null) {
                Debug.LogError("Need GravityEngine in the scene");
                return 0;
            }
            Init();
        }
        a_phy = a_scaled / gravityEngine.GetPhysicalScale();
        float massScale = gravityEngine.massScale;
        orbitPeriod = 2f * System.Math.PI * System.Math.Sqrt(a_phy * a_phy * a_phy / ((float)centerNbody.mass * massScale)); // G=1
        return orbitPeriod;
    }

    /// <summary>
    /// Called by GravityEngine to setup physics info prior to simulation
    /// Do not call this method directly. Instead ensure this object is added to the GravityEngine
    /// either by adding in a script or to the public bodies list.
    /// </summary>
    /// <param name="physicalScale">Physical scale.</param>
    public void PreEvolve(float physicalScale, float massScale) {
		if (ecc >= 1f) {
			ecc = 0.99f;
		}
		// Convert a to physical scale so the correct orbit time results w.r.t evolving bodies
		// Don't know this until factor until GravityEngine calls here. 
		a_phy = a/physicalScale;
		// base.CalculateRotation();
		// evolution relies on orbital period to determine current mean anomoly
		orbitPeriod = GetPeriod();

		// convert phase in true anomoly into mean anomoly phase
		double phase_rad = DEG_TO_RAD * phase;
		double phase_E = 2f * System.Math.Atan(System.Math.Sqrt((1-ecc)/(1+ecc))* System.Math.Tan(phase_rad/2f));
		mean_anomoly_phase = phase_E - ecc * System.Math.Sin(phase_E);
	}
		
	/// <summary>
	/// Called from the GravityEngine on FixedUpdate cycles to determine current position of body given
	/// the physics time evolution only when mode=KEPLERS_EQN.
	///
	/// This routine updates the game object position in game space and physics space. 
	///
	/// Do not call this method directly. 
	/// </summary>
	/// <param name="physicsTime">Physics time.</param>
	/// <param name="physicalScale">Physical scale.</param>
	/// <param name="r">Reference to array into which new position is placed.</param>
	public void Evolve (double physicsTime, float physicalScale, ref float[] r_new) {
		//  There is no simple expression for the position in an Elliptical orbit as a function of time. 
		//	(Expressions exist that describe position as a function of the angle in the orbit, but determining
		//	angle as a function of time results in an expression that is not elementary - Kepler's equation). 
		//			
		//	Here we follow the excellent development of the equations in "Gravity" (Poisson and Will, 2014). 
		//
		// following Gravity (Poisson & Will) Box 3.1
		// 1) First use Newton's root-finding method to determine eccentic anomoly u for a given t
		// The mean anomoly (angle from center of ellipse if motion was circular and uniform) is determined
		// from the orbit period and time evolved so far. 
		int loopCount = 0; 
		// mean_anomoly is the angle around the circle of radius a if the body moved uniformly
		// (which it does not)
		double mean_anomoly = 2f * (System.Math.PI * physicsTime/orbitPeriod);
		mean_anomoly += mean_anomoly_phase;
        double u = mean_anomoly; // seed with mean anomoly
        double u_next = 0;
		const int LOOP_LIMIT = 20;
		while(loopCount++ < LOOP_LIMIT) {
			// this should always converge in a small number of iterations - but be paranoid
			u_next = u + (mean_anomoly - (u - ecc * System.Math.Sin(u)))/(1 - ecc * System.Math.Cos(u));
			if (System.Math.Abs(u_next - u) < 1E-5)
				break;
			u = u_next;
		}
        if (loopCount >= LOOP_LIMIT) {
            u_next = u + (mean_anomoly - (u - ecc * System.Math.Sin(u))) / (1 - ecc * System.Math.Cos(u));
            Debug.LogWarning(string.Format("Failed to converge (use best estimate) err={0:E5}", System.Math.Abs(u_next - u))); 
            // keep going anyway
        }
        // 2) eccentric anomoly is angle from center of ellipse, not focus (where centerObject is). Convert
        //    to true anomoly, f - the angle measured from the focus. (see Fig 3.2 in Gravity) 
        double cos_f = (System.Math.Cos(u) - ecc)/(1 - ecc * System.Math.Cos(u));
        double sin_f = (System.Math.Sqrt(1 - ecc*ecc) * System.Math.Sin (u))/(1 - ecc * System.Math.Cos(u));
        double r = a_phy * (1f - ecc*ecc)/(1f + ecc * cos_f);
		position = new Vector3( (float)(r * cos_f), (float)(r * sin_f), 0);

        // Need a precise (right now) value for centerBody for best results, get from the engine

		// move from XY plane to the orbital plane and scale to world space
		// orbit position is WRT center
		position =  ellipse_orientation * position;
        // fill in r. NBE will use this position.
        r_new[0] = position.x;
		r_new[1] = position.y;
		r_new[2] = position.z;
		// update object position in world space
		position = physicalScale * position + GravityEngine.Instance().GetScaledPosition(centerNbody);

        // determine velocity 
        // (C&P from above for performance)
        // Murray and Dermot 
        // (2.26)
        // This should really be (M+m), but assume m << M
        // (massScale is added in at the GE level)
        double n = System.Math.Sqrt((float)(centerNbody.mass * GravityEngine.Instance().massScale) / (a_phy * a_phy * a_phy));
        // (2.36)
        double denom = System.Math.Sqrt(1f - ecc * ecc);
        double xdot = -1f * n * a_phy * sin_f / denom;
        double ydot = n * a_phy * (ecc + cos_f) / denom;
        Vector3 v_xy = new Vector3((float) xdot, (float) ydot, 0);
        velocityScaled = ellipse_orientation * v_xy + GravityEngine.Instance().GetScaledVelocity(centerObject);
        phase = NUtils.AngleFromSinCos(sin_f, cos_f) * DEG_TO_RAD;
    }


    public void Move(Vector3 moveBy) {
        position += moveBy;
    }

    public void GEUpdate(GravityEngine ge) {
        nbody.GEUpdate(position, velocityScaled, ge);
    }

    public Vector3 GetVelocity() {
        return velocityScaled;
    }

    public Vector3 GetPosition() {
        return position;
    }

	/// <summary>
	/// Convert Mean Anomoly to True Anomoly for an ellipse with eccentricity e. 
	/// </summary>
	/// <returns>True Anomoly in degrees.</returns>
	/// <param name="m">Mean Anomoly. (degrees)</param>
	/// <param name="e">Eccentricty.</param>
	public static double MeanToTrueAnomoly(double m, double e) {
		int loopCount = 0;
        double u = m * DEG_TO_RAD; // seed with mean anomoly
        double u_next = 0;
		// some extreme comet orbits (e.g. Halley) need a LOT of iterations
		const int LOOP_LIMIT = 200;
		while(loopCount++ < LOOP_LIMIT) {
			// this should always converge in a small number of iterations - but be paranoid
			u_next = u + (m - (u - e * System.Math.Sin(u)))/(1 - e * System.Math.Cos(u));
			if (System.Math.Abs(u_next - u) < 1E-6)
				break;
			u = u_next;
		}
		if (loopCount >= LOOP_LIMIT)	
			Debug.LogError("Failed to converge u_n=" + u_next);	// keep going anyway

		// 2) eccentric anomoly is angle from center of ellipse, not focus (where centerObject is). Convert
		//    to true anomoly, f - the angle measured from the focus. (see Fig 3.2 in Gravity) 
		double cos_f = (System.Math.Cos(u) - e)/(1f - e * System.Math.Cos(u));
        double sin_f = (System.Math.Sqrt(1 - e*e) * System.Math.Sin (u))/(1f - e * System.Math.Cos(u));
        double f_deg = NUtils.AngleFromSinCos(sin_f, cos_f) * DEG_TO_RAD;
        //Debug.Log("m=" + m + " E=" + u * System.Math.Deg2Rad + " f=" + f_deg);
        return f_deg;
	}

    // Allow either P or A as params. Only one can be specified in editor at a time since they are related.
    // If the editor changes the value, the keep the related variable in sync.
    private void UpdateOrbitParams() {
        if (paramBy == ParamBy.AXIS_A) {
            p = a * (1 - ecc);
        } else if (paramBy == ParamBy.CLOSEST_P) {
            a = p / (1 - ecc);
        }
    }

    /// <summary>
    /// Apply scale to the orbit. This is used by the inspector scripts during
    /// scene setup. Do not use at run-time.
    /// </summary>
    /// <param name="scale">Scale.</param>
    public void ApplyScale(float scale) {
		if (paramBy == ParamBy.AXIS_A){
			a_scaled = a * scale;
		} else if (paramBy == ParamBy.CLOSEST_P) {
			p_scaled = p * scale; 
			a_scaled = p * scale/(1-ecc);
		}
		UpdateOrbitParams();
        nbody = GetComponent<NBody>();
        SetInitialPosition(nbody);
	}

	public void Log(string prefix) {
		Debug.Log(string.Format("orbitEllipse: {0} a={1} e={2} i={3} Omega={4} omega={5} phase={6}", 
								prefix, a, ecc, inclination, omega_uc, omega_lc, phase));
	}
}
