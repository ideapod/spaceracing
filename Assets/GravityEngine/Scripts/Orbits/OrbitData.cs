using UnityEngine;
using System.Collections;

/// <summary>
/// Orbit data.
/// Hold the traditional orbit parameters for an elliptic/hyperbolic orbit. 
///
/// Provide utility methods to derive the orbit parameters from the position, velocity and centerBody
/// of an NBody object pair. This orbit prediction is based only on the two bodies (and assumes the central
/// body mass dominates) - the influence of other bodies in the scene may result in an actual path that is
/// different from the predicted path. 
/// </summary>
[System.Serializable]
public class OrbitData  {

    public const float OD_ZERO_LIMIT = 0.001f;

	// Orbit parameters (user control via FixedEllipseEditor)
	// These parameters are in world space. 
	//! eccentricity (0..1, 0=circle, 1=linear)
	public float ecc; 			

	// Allow EITHER a or p to specify size (JPL data uses a for planets and p for comets, so both are useful)
	//! semi-major axis - based on paramBy user can specify a OR p. a = p/(1-ecc)
	public float a = 10f; 			
	//! perihelion 
	public float perihelion; 			
	//! "longitude of ascending node" - angle from x-axis to line from focus to pericenter
	public float omega_uc; 		
	//! "argument of perienter" - angle from ascending node to pericenter
	public float omega_lc; 		
	//! inclination (degrees!)
	public float inclination; 	

	//! initial TRUE anomoly (angle wrt line from focus to closest approach)
	public float phase;

    //! period of orbit in game seconds (use GravityScalar.GameSecToWorldTime for world units)
    public float period;

 	//! time to periapsis (point of closest approach) in game seconds (use GravityScalar.GameSecToWorldTime for world units)
	public float tau; 	

	//! Hyperbola - initial distance from focus
	public float r_initial = 10f; 	

	public NBody centralMass;

	public NBody nbody;

	//! Scaled mass of central body
	public float mu;

    private const float MIN_ANGULAR_MOMTM = 0.001f;

    private const float MIN_VECTOR_LEN = 0.01f;

    // Empty constructor
    public OrbitData() {
        // empty
    }

    /// <summary>
    /// Construct an orbit data routine from an existing orbit ellipse by copying the orbital elements
    /// </summary>
    /// <param name="orbitEllipse"></param>
    public OrbitData(OrbitEllipse orbitEllipse) {
        a = orbitEllipse.a_scaled;
        omega_lc = orbitEllipse.omega_lc;
        omega_uc = orbitEllipse.omega_uc;
        inclination = orbitEllipse.inclination;
        ecc = orbitEllipse.ecc;
        phase = orbitEllipse.phase;
        centralMass = orbitEllipse.centerObject.GetComponent<NBody>();
        mu = centralMass.mass * GravityEngine.Instance().massScale;
    }

    /// <summary>
    /// Computes the orbit parameters for a specified velocity with respect to a central body.
    /// This method assumes the orbit is not a Kepler-mode orbit ellipse. If it is (or not sure) then use
    /// SetOrbit() below.
    /// </summary>
    /// <param name="forNbody"></param>
    /// <param name="aroundNBody"></param>
    /// <param name="velocity"></param>
    public void SetOrbitForVelocity(NBody forNbody, NBody aroundNBody, Vector3 velocity) {
        nbody = forNbody;
        centralMass = aroundNBody;
        Vector3 v_rel = velocity - aroundNBody.vel_phys;
        float velSq = v_rel.sqrMagnitude;
        Vector3 r_for = GravityEngine.Instance().GetPhysicsPosition(forNbody);
        Vector3 r_around = GravityEngine.Instance().GetPhysicsPosition(aroundNBody);
        Vector3 r_vec = r_for - r_around;
        float r = Vector3.Magnitude(r_vec);
        mu = aroundNBody.mass * GravityEngine.Instance().massScale;

        // Vallado Algorithm 9 p113
        Vector3 ecc_vec = (1 / mu) * ((velSq - mu / r) * r_vec - Vector3.Dot(r_vec, v_rel) * v_rel);
        ecc = ecc_vec.magnitude;

        // vis-visa for ellipse (Roy 4.36) [compare 4.90 for hyperbola]
        if (ecc < 1.0f) {
            EllipseOrbitForVelocity(forNbody, aroundNBody, r_vec, v_rel);
        } else {
            HyperOrbitForVelocity(forNbody, aroundNBody, r_vec, v_rel);
        }
        //Debug.Log(string.Format("SetOrbit: a={0} perih={1} e={2} i={3} Omega={4} omega={5} r_initial={6}",
        //                a, perihelion, ecc, inclination, omega_uc, omega_lc, r_initial));
    }


    /// <summary>
    /// Computes the orbit parameters for a specified velocity (from NBody) with respect to a central body.
    /// This method assumes the orbit is not a Kepler-mode orbit ellipse. If it is (or not sure) then use
    /// SetOrbit() below.
    /// </summary>
    /// <param name="forNbody"></param>
    /// <param name="aroundNBody"></param>
    public void SetOrbitForVelocity(NBody forNbody, NBody aroundNBody) {
		nbody = forNbody;
        nbody.UpdateVelocity();
        centralMass = aroundNBody;
        aroundNBody.UpdateVelocity();
        SetOrbitForVelocity( forNbody, aroundNBody, forNbody.vel_phys);
    }

    /// <summary>
    /// SetOrbit
    /// Determine orbit params from an attached orbit component (if Kepler) otherwise use the velocity to 
    /// determine the orbit
    /// </summary>
    /// <param name="forNbody"></param>
    /// <param name="aroundNBody"></param>
    public void SetOrbit(NBody forNbody, NBody aroundNBody) {
        nbody = forNbody;
        centralMass = aroundNBody;
        // is this a Kepler body
        if (nbody.engineRef.fixedBody != null) {
            OrbitEllipse orbitEllipse = nbody.GetComponent<OrbitEllipse>();
            if (orbitEllipse != null) {
                ecc = orbitEllipse.ecc;
                a = orbitEllipse.a * GravityEngine.Instance().GetLengthScale();
                omega_lc = orbitEllipse.omega_lc;
                omega_uc = orbitEllipse.omega_uc;
                inclination = orbitEllipse.inclination;
                period = CalcPeriod(a, aroundNBody);
                // TODO:  tau
                phase = orbitEllipse.phase;
                return;
            }
            OrbitHyper orbitHyper = nbody.GetComponent<OrbitHyper>();
            if (orbitHyper != null) {
                ecc = orbitHyper.ecc;
                perihelion = orbitHyper.perihelion * GravityEngine.Instance().GetLengthScale();
                omega_lc = orbitHyper.omega_lc;
                omega_uc = orbitHyper.omega_uc;
                inclination = orbitHyper.inclination;
                // need phase, tau, period
                return;
            }
        }
        SetOrbitForVelocity(forNbody, aroundNBody);
    }

    /// <summary>
    /// Compute the axis of the orbit by taking r x v and normalizing. 
    /// </summary>
    /// <returns></returns>
    public static Vector3d GetAxis(NBody forNbody, NBody aroundNBody) {
        GravityEngine ge = GravityEngine.Instance();
        Vector3d r = ge.GetPositionDoubleV3(forNbody) - ge.GetPositionDoubleV3(aroundNBody);
        Vector3d v = GravityEngine.Instance().GetVelocityDoubleV3(forNbody) - ge.GetVelocityDoubleV3(aroundNBody);
        return Vector3d.Normalize( Vector3d.Cross(r, v));
    }

    private float CalcPeriod(float a, NBody centerMass) {
        // this equation has a G but we set G=1
        float _period = 2f * Mathf.PI * Mathf.Sqrt(a * a * a) / Mathf.Sqrt(centerMass.mass);
        return GravityScaler.ScaleToGameSeconds(_period); 
    }

    private void EllipseOrbitForVelocity(NBody forNbody, NBody aroundNBody, Vector3 r_vec, Vector3 velocity) {
        // Chobotov 3rd Ed, Sec 4.9 p66
        Vector3 vel = velocity;
        float velSq = vel.sqrMagnitude;

        float r = Vector3.Magnitude(r_vec);
        float mu = aroundNBody.mass * GravityEngine.Instance().massScale;

        // a - semi-major axis
        a = 1f / (2f / r - velSq / mu);

        // Debug.Log("r=" + r + " r_vec=" + r_vec + " vel= " + vel + " forVel=" + forNbody.vel_phys + " aroundVel=" + aroundNBody.vel_phys);

        // W is the normal to orbital plane. The angle to z-axis determines inclination. 
        Vector3 W = Vector3.Normalize(Vector3.Cross(r_vec, vel));   // 4.118
        float inclination_rad = Mathf.Acos(Vector3.Dot(W, Vector3.forward)); // forward = z-axis
        inclination = Mathf.Rad2Deg * inclination_rad;

        // Debug.LogFormat("incl={0} W={1}", inclination, W);
        // TODO: Pass in ecc instead of re-compute
        Vector3 ecc_vec = 1 / mu * ((velSq - mu / r) * r_vec - Vector3.Dot(r_vec, vel) * vel); // 4.120
        ecc = Vector3.Magnitude(ecc_vec);

        omega_lc = 0f;
        omega_uc = 0f;
        // Omega is angle to line of nodes and applies when incl != 0
        // omega is angle to pericenter and is zero if ecc = 0

        // Omega 4.122 (i-hat is x-axis)
        // N is the unit vector indicating line of nodes. It is only defined if there is non-zero inclination
        // Use Sin to pick up 0 and 180. 
        Vector3 N = Vector3.zero;
        if (Mathf.Abs(Mathf.Sin(inclination_rad)) > 1E-4) {
            // inclination 
            N = Vector3.Normalize(Vector3.Cross(Vector3.forward, W));
            float cos_Omega = Vector3.Dot(Vector3.right, N);
            float sin_Omega = Vector3.Dot(Vector3.Cross(Vector3.right, N), Vector3.forward);
            omega_uc = Mathf.Rad2Deg * NUtils.AngleFromSinCos(sin_Omega, cos_Omega);

            // omega
            if (ecc > MIN_VECTOR_LEN) {
                float cos_omega = Vector3.Dot(N, ecc_vec) / ecc;
                float sin_omega = Vector3.Dot(Vector3.Cross(N, ecc_vec), W) / ecc;
                omega_lc = Mathf.Rad2Deg * NUtils.AngleFromSinCos(sin_omega, cos_omega); ;
                // Debug.LogFormat(" sin={0} cos={1} omega={2} e.x={3} e.y={4}", sin_omega, cos_omega, omega_lc, ecc_vec.x, ecc_vec.y);
            }
        } else {
            // inclination = 0 or 180. No line of nodes, no Omega. 
            if (ecc > MIN_VECTOR_LEN) {
                float sin_omega = ecc_vec.y / ecc;
                float cos_omega = ecc_vec.x / ecc;
                omega_lc = Mathf.Rad2Deg * NUtils.AngleFromSinCos(sin_omega, cos_omega);
                // Emperical HACK - is this really an orbit renderer bug?
                if (inclination > 90) {
                    omega_lc = -omega_lc;
                }
                // Debug.LogFormat(" sin={0} cos={1} omega={2} e.x={3} e.y={4}", sin_omega, cos_omega, omega_lc, ecc_vec.x, ecc_vec.y);
            }
        }
        // prefer 0 to 360
        if (Mathf.Abs(omega_lc-360f) < 1E-4) {
            omega_lc = 0f; 
        }

        // phase
        // Need a special case for e=0
        float sin_theta = 0f;
        float cos_theta = 0f;
        float er = ecc * r;
        // ansatz: if ecc is zero need to find phase wrt N (if inclination), else wrt x-axis
        if (ecc < MIN_VECTOR_LEN) {
            if (Mathf.Abs(inclination) > MIN_VECTOR_LEN) {
                ecc_vec = N;
            } else {
                ecc_vec = Vector3.right;
            }
            er = r;
        }
        cos_theta = Vector3.Dot(ecc_vec, r_vec) / er;
        sin_theta = Vector3.Dot(Vector3.Cross(ecc_vec, r_vec), W) / er;
        phase = Mathf.Rad2Deg * NUtils.AngleFromSinCos(sin_theta, cos_theta); ;

        // Determine time to periapsis (M&D 2.140)
        // Eqn assumes E=0 when t=tau
        float E = Mathf.Acos((1f - r / a) / ecc);
        // this equation has a G but we set G=1
        period =  CalcPeriod(a, aroundNBody);   // calc will also scale

        float M = (E - ecc * Mathf.Sin(E));
        tau = M * Mathf.Sqrt(a * a * a) / Mathf.Sqrt(aroundNBody.mass);
        tau = GravityScaler.ScaleToGameSeconds(tau);
        // tau is giving time to/from apoapsis, need to find out which
        float vdotr = Vector3.Dot(vel, r_vec);
        if (vdotr > 0) {
            tau = period - tau;
        }
    }

    private void HyperOrbitForVelocity(NBody forNbody, NBody aroundNBody, Vector3 r_vec, Vector3 velocity) {
        // Murray and Dermott Ch2.8 (2.134) - (2.140)
        Vector3 vel = velocity - aroundNBody.vel_phys;
        float velSq = vel.sqrMagnitude;

		float r = Vector3.Magnitude(r_vec);
		r_initial = r;
		float mu = aroundNBody.mass * GravityEngine.Instance().massScale;
        // Vallado p114
        float xi = velSq / 2 - mu / r;
        a = 0.5f * mu / xi;

		// Determine angular momentum, h
		Vector3 h_vec = Vector3.Cross(r_vec, vel);
		float h = Vector3.Magnitude(h_vec);
		//Debug.Log("h_vec = " + h_vec + " r_vec=" + r_vec + " v=" + forNbody.vel + " planet pos=" + forNbody.transform.position);

		perihelion = a*(ecc-1);

        // p = a*(1-ecc^2);

		// inclination (0..180 so Acos is unambiguous) (2.136)
		float inclRad = Mathf.Acos(h_vec.z/h);

		// Omega_uc - only relevant if there is non-zero inclination
		float omega_ucRad = 0f;
		float sinOmega = 0f;
		float cosOmega = 1f;
		bool inclNoneZero = Mathf.Abs(Mathf.Sin(inclRad)) > 1E-5;
		if (inclNoneZero) {
			sinOmega = h_vec.x/(h*Mathf.Sin(inclRad));
			cosOmega = -h_vec.y/(h*Mathf.Sin(inclRad));
		} else if (h_vec.z < 0) {
			// if incl = 180
			sinOmega *= -1f;
			cosOmega *= -1f;
		}
		omega_ucRad = NUtils.AngleFromSinCos( sinOmega, cosOmega);
//		Debug.Log("sinOmega=" + sinOmega + " cosOmega=" + cosOmega + " omega_ucRad=" + omega_ucRad + " inclRad=" + inclRad
//			+ " inclNonZero=" + inclNoneZero + " h_vec=" + h_vec);

		float f = 0; 
		if (ecc > 1E-3) {
			// Like M&D 2.31 but for e^2-1 not 1-e^2
			float rdot = Mathf.Sqrt(velSq - h*h/(r*r));
			if (Vector3.Dot(r_vec, forNbody.vel_phys) < 0) {
				rdot *= -1f;
			} 
			float cos_f = (1/ecc)*(a*(ecc*ecc-1)/r - 1f);
			float sin_f = rdot*a*(ecc*ecc-1f)/(h*ecc);
			f = NUtils.AngleFromSinCos(sin_f, cos_f);
		} 

		// (2.6.17)
		// sin(w+f) = Z/(r sin(i))
		// cos(w+f) = [X cos(Omega) + Y sin(Omega)]/r
		// when sin(i) = 0, will have Z=0 and w + f = 0
		float sin_of = 0f;
		float cos_of = 1f;
		float omega_lcRad = 0f;
		if (inclNoneZero) {
			sin_of = r_vec.z/(r*Mathf.Sin(inclRad));
			cos_of = (r_vec.x * cosOmega + r_vec.y * sinOmega)/r;
			omega_lcRad = (NUtils.AngleFromSinCos(sin_of, cos_of) - f);
		} else {
			float sin_theta = r_vec.y/r;
			float cos_theta = r_vec.x/r;
			float theta = NUtils.AngleFromSinCos(sin_theta, cos_theta);
//			Debug.Log("r_vec=" + r_vec + " theta=" + theta + " f=" + f + " (theta - f)=" + (theta-f));
			if (inclRad < Mathf.PI/2f) {
				omega_lcRad = theta - f; 
			} else {
				if (h_vec.z > 0) {
					omega_ucRad = f - theta; 
				} else {
					// Rotated by 180 around x, so flip sign of y
					sin_theta = -r_vec.y/r;
					cos_theta = r_vec.x/r;
					theta = NUtils.AngleFromSinCos(sin_theta, cos_theta);
					omega_ucRad = f - theta; 
				}
			}
		}
		// Ellipse wants degrees
		omega_lc = Mathf.Rad2Deg * omega_lcRad;
		if (omega_lc < 0)
			omega_lc += 360f;
		omega_uc = Mathf.Rad2Deg * omega_ucRad;
		if (omega_uc < 0)
			omega_uc += 360f;
		phase = Mathf.Rad2Deg * f; 
		inclination = Mathf.Rad2Deg * inclRad;
//		Debug.Log("cosof=" + cos_of + " sinof=" + sin_of + " of=" + NUtils.AngleFromSinCos(sin_of, cos_of) + " f=" + f);

	}

    public string LogString() {
        return string.Format("a={0} ecc={1} incl={2} Omega={3} omega={4} phase ={5} period={6} tau={7}", 
                        a, ecc, inclination, omega_uc, omega_lc, phase, period, tau);
    }

    /// <summary>
    /// Get the physics position for a specified phase for the ellipse defined by this OrbitData.
    /// </summary>
    /// <param name="phaseDeg"></param>
    /// <returns></returns>
    public Vector3 GetPhysicsPositionforEllipse(float phaseDeg) {

        // C&P from EllipseBase - make common someday
        Quaternion ellipse_orientation = Quaternion.AngleAxis(omega_uc, GEConst.zunit) *
                      Quaternion.AngleAxis(inclination, GEConst.xunit) *
                      Quaternion.AngleAxis(omega_lc, GEConst.zunit);

        float phaseRad = phaseDeg * Mathf.Deg2Rad;
        // position object using true anomoly (angle from  focus)
        // a is really a_scaled when constructed from OrbitEllipse (so scaling has been done)
        float r = a * (1f - ecc * ecc) / (1f + ecc * Mathf.Cos(phaseRad));

        Vector3 pos = new Vector3(r * Mathf.Cos(phaseRad), r * Mathf.Sin(phaseRad), 0);
        // move from XY plane to the orbital plane
        Vector3 new_p = ellipse_orientation * pos;
        // orbit position is WRT center. Could be adding dynamically to an object in motion, so need current position. 
        Vector3 centerPos = Vector3.zero;
        // used by widgets - so need to get explcitly
        if (centralMass.engineRef != null) {
            centerPos = GravityEngine.Instance().GetPhysicsPosition(centralMass);
        } else {
            // setup - not yet added to GE
            centerPos = centralMass.initialPhysPosition;
        }
        return new_p + centerPos;
    }

    /// <summary>
    /// Determine velocity in physics units at the indicated phase angle (in degrees)
    /// 
    /// </summary>
    /// <param name="phaseDeg"></param>
    /// <returns></returns>
    public Vector3 GetPhysicsVelocityForEllipse(float phaseDeg) {
        /// Uses Vallado, Algorithm 10 for (x,y) plane and then rotates into place

        Quaternion ellipse_orientation = Quaternion.AngleAxis(omega_uc, GEConst.zunit) *
                       Quaternion.AngleAxis(inclination, GEConst.xunit) *
                       Quaternion.AngleAxis(omega_lc, GEConst.zunit);
        float p = a * (1f - ecc * ecc);
        float phaseRad = phaseDeg * Mathf.Deg2Rad;
        float vx = -Mathf.Sqrt(mu / p) * Mathf.Sin(phaseRad);
        float vy = Mathf.Sqrt(mu / p) * (ecc + Mathf.Cos(phaseRad));
        return ellipse_orientation * new Vector3(vx, vy, 0);
    }

}
