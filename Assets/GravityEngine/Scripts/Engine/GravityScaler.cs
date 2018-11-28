using UnityEngine;
using System.Collections;

public class GravityScaler : MonoBehaviour  {

	/// <summary>
	/// Units.
	/// Gravity N-Body simulations use units in which G=1 (G, the gravitational constant). 
	/// 
	/// From unit analysis of F = m_1 a = (G m_1 m_2)/r^2 we get: T = (L^3/G M)^(1/2) with
	///  T = time, L = length, M = Mass, G=Newton's constant = 6.67E-11.
	///
	///  SI (m,kg) => T=1.224E5 sec
	///  Solar (AU, 1E24 kg) => T = 7.083E9 sec.
	///
	/// To control game time, mass is rescaled in the physics engine to acheive the desired result. 
	/// Initial velocities are also adjusted to appropriate scale. 
	/// </summary>
	public const float G = 6.67408E-11f;  // m^3 kg^(-1) sec^(-2)

	public enum Units { DIMENSIONLESS, SI, ORBITAL, SOLAR }; //, STELLAR};
	public Units units;
	private static string[] lengthUnits =   {"DL", "m", "km", "AU", "light-year"};
	private static string[] massUnits =     {"DL", "kg", "1E24 kg", "1E24 kg", "Msolar"};
	private static string[] velocityUnits = {"DL", "m/s", "km/hr", "km/s", "km/s"};

	private static float velocityScale;

    public const float M_PER_KM = 1000f;
    public const float M_PER_AU = 1.49598E+11f;
    public const float SEC_PER_YEAR = 3600f*24f*365.25f;
    public const float KM_HOUR_TO_M_SEC = M_PER_KM/SEC_PER_HOUR;

    public const float KM_SEC_TO_AU_YR = 0.210949527f;

	// Use T = Sqrt(L^3/(G M)) with all units and G in m/s/kg.
	// This gives the time unit that results from the choice of length/mass and the imposition of G=1
	// in the integrators. 
	// m/kg/sec
	private const float G1_TIME_SI = 122406.4481f;
	// km/ 1E24kg/ hr
	public const float SEC_PER_HOUR = 3600f;
	private const float G1_TIME_ORBIT = 0.003870832f;

	// SOLAR Units
	private const float G1_TIME_SOLAR = 7082595090f;

    public static float game_sec_per_phys_sec = 1f;

    // i.e. time is in units of approx. 3 ms
    // If we want Unity time in game sec. per physics hour then we need to 
    // timescale = game sec. per physics hour

    /// <summary>
    /// Determine and set the required massScale given:
    /// - inspector provided timeScale
    /// - inspector provided lengthScale
    /// - choice of units in inspector (Dimensionless, SI, Orbital, Solar)
    /// 
    /// The integrators in GE all assume G=1. 
    /// The time evolution in GE is controlled by scaling mass to acheive the desired speed. 
    /// 
    /// </summary>
    /// <param name="units"></param>
    /// <param name="timeScale"></param>
    /// <param name="lengthScale"></param>
    public static void UpdateTimeScale(Units units, float timeScale, float lengthScale) {

		// time unit size (sec.) for G=1 and units given
		float time_g1 = 1f; 
		// Number of physics seconds per Unity second, given the units and timeScale
		// Convert all cases to SI units - this is the units for G
		// Adjust the SI units by the Unity scale factors. 
		switch(units) {
		case Units.DIMENSIONLESS:
			// mass scale is controlled via inspector for dimensionless units
			time_g1 = G1_TIME_SI;
			game_sec_per_phys_sec = 1f/timeScale;
			velocityScale = 1f;
			#pragma warning disable 162		// disable unreachable code warning
			if (GravityEngine.DEBUG) {
				Debug.Log("SetMassScale: Time G1 = " + time_g1 + 
						" game_sec_per_phys_sec=" + game_sec_per_phys_sec + 
						" velScale=" + velocityScale );
			}
			#pragma warning restore 162
			return;
		case Units.SI:
			time_g1 = G1_TIME_SI;
			game_sec_per_phys_sec = 1f/timeScale;
			velocityScale = lengthScale/timeScale;
			break;
		case Units.ORBITAL:
			time_g1 = G1_TIME_ORBIT;
			game_sec_per_phys_sec = SEC_PER_HOUR/timeScale;
			velocityScale = lengthScale/timeScale;
			break;
		case Units.SOLAR:
			time_g1 = G1_TIME_SOLAR;
			game_sec_per_phys_sec = SEC_PER_YEAR/timeScale;
			velocityScale = KM_SEC_TO_AU_YR * lengthScale/timeScale;
			break;
		default:
			Debug.LogError("Unsupported units");
			break;
		}

		// The length scale chosen for the scene is now applied.
		// Convert to the designated length scale (convert from m to Unity units)
		float time_unity = time_g1 / Mathf.Sqrt(lengthScale*lengthScale*lengthScale);

		// timeScale indicates game seconds per physics sec. 

		// The masses do not affect position in scene, so instead of adjusting raw masses of all NBody objects this
		// adjusment is done to the physics copy in GE as they are added. See GravityEngine:SetupOneGameObject() [private]
		float mu = game_sec_per_phys_sec*game_sec_per_phys_sec/(time_unity*time_unity);
		GravityEngine.Instance().massScale = mu;

		#pragma warning disable 162		// disable unreachable code warning
		if (GravityEngine.DEBUG) {
			Debug.Log("SetMassScale: Time G1 = " + time_g1 + 
					" time_unity=" + time_unity +
                    " timeScale=" + timeScale +
					" game_sec_per_phys_sec=" + game_sec_per_phys_sec + 
					" massScale=" + mu + " velScale=" + velocityScale );
		}
		#pragma warning restore 162

	}

    public static float GetGameSecondPerPhysicsSecond()
    {
        return game_sec_per_phys_sec;
    }

    [System.Obsolete("Use the better named ScaleToGameSeconds")]
    public static float ScalePeriod(float period)
    {
        // do not need to scale length - already happened in determination of "a"
        float massScale = GravityEngine.Instance().massScale;
        return period / Mathf.Sqrt(massScale) ;
    }

    /// <summary>
    /// Scale a physics time to Unity game seconds
    /// </summary>
    /// <param name="time"></param>
    /// <returns></returns>
    public static float ScaleToGameSeconds(float time) {
        float massScale = GravityEngine.Instance().massScale;
        return time / Mathf.Sqrt(massScale);
    }

    public static float GameSecToWorldTime(float gameSec) {
        return gameSec * game_sec_per_phys_sec;
    }

    /// <summary>
    /// Modify the physics engine acceleration to appear in world scale units. 
    /// </summary>
    /// <param name="a"></param>
    /// <param name="lengthScale"></param>
    /// <param name="timeScale"></param>
    /// <returns></returns>
    public static Vector3 ScaleAcceleration(Vector3 a, float lengthScale, float timeScale)
    {
        return  (timeScale * timeScale) * a / lengthScale ;
    }

	public static float GetVelocityScale() {
		return velocityScale;
	}

	/// <summary>
	/// Changes the length scale of all NBody objects in the scene due to a change in the inspector.
	/// Find all NBody containing objects.
	/// - independent objects are rescaled
	/// - orbit based objects have their primary dimension adjusted (e.g. for ellipse, a)
	///   (these objects are scalable and are asked to rescale themselves)
	///
	/// Not intended for run-time use.
	/// </summary>
	public static void ScaleScene(Units units, float lengthScale) {

		// ensure velocity scale is determined
		// find everything with an NBody. Rescale will ensure only independent NBodies are rescaled. 
		NBody[] nbodies = FindObjectsOfType<NBody>();
		foreach (NBody nbody in nbodies) {
			ScaleNBody(nbody, units, lengthScale);
		}
	}

	/// <summary>
	/// Scales the N body using the provided units and length scale.
	/// </summary>
	/// <param name="nbody">Nbody.</param>
	/// <param name="units">Units.</param>
	/// <param name="timeScale">Time scale.</param>
	/// <param name="lengthScale">Length scale.</param>
	public static void ScaleNBody(NBody nbody, Units units, float lengthScale) {
		// If there is an IOrbitScaler - use it instead
		IOrbitScalable iorbit = nbody.GetComponent<IOrbitScalable>();
		if (iorbit != null) {
			iorbit.ApplyScale(lengthScale);
		} else {
			if (units == GravityScaler.Units.DIMENSIONLESS && lengthScale == 1f) {
				// Backwards compatibility with pre 1.3 GE
				nbody.initialPos = nbody.transform.position;
			}
			nbody.ApplyScale(lengthScale, velocityScale);	
		}
	}

    //-----------CONVERSION HELPERS-------------------

    /// <summary>
    /// Provide conversion factor to convert acceleration in SI units to the
    /// GE units specified in the argument. 
    /// 
    /// This routine references the currently set units in the GE
    /// </summary>
    /// <returns>Conversion Factor</returns>
    public static double AccelSItoGEUnits()
    { 
        double conversion = 1.0;
        switch(GravityEngine.Instance().units)
        {
            case Units.ORBITAL:
                conversion = SEC_PER_HOUR * SEC_PER_HOUR / M_PER_KM;
                break;

            case Units.SOLAR:
                Debug.LogError("Implement me");
                break;

            case Units.DIMENSIONLESS:
            case Units.SI:
                // nothing to do 
                break;
        }
        return conversion;
    }

    /// <summary>
    /// Determine the factor to convert acceleration to SI units. 
    /// </summary>
    /// <returns></returns>
    public static double AccelGEtoSIUnits()
    {
        double conversion = 1.0;
        switch (GravityEngine.Instance().units)
        {
            case Units.ORBITAL:
                conversion = M_PER_KM /(SEC_PER_HOUR * SEC_PER_HOUR) ;
                break;

            case Units.SOLAR:
                Debug.LogError("Implement me");
                break;

            case Units.DIMENSIONLESS:
            case Units.SI:
                // nothing to do 
                break;
        }
        return conversion;
    }

    /// <summary>
    /// Determine the factor to convert scaled position to SI units. 
    /// </summary>
    /// <returns></returns>
    public static double PositionScaletoSIUnits() {
        double conversion = 1.0;
        switch (GravityEngine.Instance().units) {
            case Units.ORBITAL:
                // KM to M
                conversion = M_PER_KM;
                break;

            case Units.SOLAR:
                // AU to m
                conversion = M_PER_AU;
                break;

            case Units.DIMENSIONLESS:
            case Units.SI:
                // nothing to do 
                break;
        }
        return conversion;
    }

    /// <summary>
    /// Determine the factor to convert scaled velocity to SI units. 
    /// </summary>
    /// <returns></returns>
    public static double VelocityScaletoSIUnits() {
        double conversion = 1.0;
        switch (GravityEngine.Instance().units) {
            case Units.ORBITAL:
                // km/hr to m/s
                conversion = M_PER_KM / SEC_PER_HOUR;
                break;

            case Units.SOLAR:
                // km/sec to m/s
                conversion = M_PER_KM;
                break;

            case Units.DIMENSIONLESS:
            case Units.SI:
                // nothing to do 
                break;
        }
        return conversion;
    }

    /// <summary>
    /// Provide the conversion factor to convert acceleration from GE units
    /// (e.g. km/hr^2 if ORBITAL) into the internal units. The internal units
    /// are determined based on the timeScale and lengthScale chosen. 
    /// 
    /// This routine references the currently set units in the GE
    /// </summary>
    /// <returns>Conversion factor</returns>
    public static double AccelGEtoInternalUnits()
    {
        GravityEngine ge = GravityEngine.Instance();
        return ge.lengthScale / (ge.timeScale * ge.timeScale);
    }

    //-------------STRING HELPERS----------------------

	/// <summary>
	/// Return the string indicating the length units in use by the gravity engine. 
	/// </summary>
	/// <returns>The units.</returns>
	public static string LengthUnits(Units units) {
		return lengthUnits[(int) units];
	}

	/// <summary>
	/// Return the string indicating the velocity units in use by the gravity engine. 
	/// </summary>
	/// <returns>The units.</returns>
	public static string VelocityUnits(Units units) {
		return velocityUnits[(int) units];
	}

	/// <summary>
	/// Return the string indicating the mass units in use by the gravity engine.
	/// </summary>
	/// <returns>The units.</returns>
	public static string MassUnits(Units units) {
		return massUnits[(int) units];
	}

    /// <summary>
    /// Return the world time in the selected units as a string
    /// </summary>
    /// <returns></returns>
    public static string GetWorldTime(double physTime, Units units) {
        string s;
        double time = physTime * game_sec_per_phys_sec;
        switch (units) {
            case GravityScaler.Units.DIMENSIONLESS:
                s = string.Format("{0:0000.0}", time);
                break;
            case GravityScaler.Units.ORBITAL:
                // convert to HHh MMm SSs (24 hour clock)
                // use ticks constructor (100ns intervals is an integer tick)
                long ticTime = (long)(time * 1E7);
                System.DateTime dateTime = new System.DateTime(ticTime);
                s = dateTime.ToString("HH:mm:ss"); // 24h format
                break;
            case GravityScaler.Units.SI:
                s = string.Format("{0:0000.0}", time);
                break;
            case GravityScaler.Units.SOLAR:
                // convert to HHh MMm SSs (24 hour clock)
                // use ticks constructor (100ns intervals is an integer tick)
                long ticTimeSolar = (long)(time * 1E7);
                System.DateTime dateTimeSolar = new System.DateTime(ticTimeSolar);
                s = dateTimeSolar.ToString("yyyy:MM:dd:HH:mm"); // 24h format
                break;
            default:
                s = "unknown units";
                break;
        }
        return s;
    }

    }
