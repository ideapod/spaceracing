using UnityEngine;
using System.Collections;

public class OrbitUtils  {

	/// <summary>
	/// Calculates the Hill Radius (radius at which the secondary's gravity becomes dominant, when the 
	/// secondary is in orbit around the primary). 
	/// </summary>
	/// <returns>The radius.</returns>
	/// <param name="primary">Primary.</param>
	/// <param name="secondary">Secondary. In orbit around primary</param>
	static public float HillRadius(GameObject primary, GameObject secondary) {

		NBody primaryBody = primary.GetComponent<NBody>(); 
		NBody secondaryBody = secondary.GetComponent<NBody>(); 
		EllipseBase orbit = secondary.GetComponent<EllipseBase>();
		if ((primaryBody == null) || (secondaryBody == null) || (orbit == null)) {
			return 0;
		}
		float denom = 3f*(secondaryBody.mass + primaryBody.mass);
		if (Mathf.Abs(denom) < 1E-6) {
			return 0;
		}
		return Mathf.Pow(secondaryBody.mass/denom, 1/3f) * orbit.a_scaled * (1-orbit.ecc);

	}

    /// <summary>
    /// Get the center
    /// </summary>
    /// <param name="objectInOrbit"></param>
    /// <returns></returns>
    static public NBody GetCenterNbody(Transform objectInOrbit, GameObject centerObject) {
        // If parent has an Nbody assume it is the center
        NBody centerNbody = null;
        if (centerObject == null) {
            if (objectInOrbit.parent != null) {
                centerNbody = objectInOrbit.parent.gameObject.GetComponent<NBody>();
                if (centerNbody != null) {
                    centerObject = objectInOrbit.parent.gameObject;
                } else {
                    Debug.LogError("Parent object must have NBody attached");
                    return null;
                }
            } else {
                Debug.Log("Warning - Require a parent object (with NBody)");
                // This path when init-ed via Instantiate() script will need to 
                // call Init() explicily once orbit params and center are set
                return null;
            }
        } else {
            centerNbody = centerObject.GetComponent<NBody>();
            if (centerNbody == null) {
                Debug.LogError("CenterObject must have an NBody attached");
            }
        }
        return centerNbody;
    }
}
