using UnityEngine;
using System.Collections;

/// <summary>
/// Fixed object.
///
/// Object does not move (but it's gravity will affect others). 
///
/// Good choice for e.g. central star in a system
/// </summary>
public class FixedObject : MonoBehaviour, IFixedOrbit {

    private NBody nbody;

    private Vector3 phyPosition; 

    public void Start() {
        nbody = GetComponent<NBody>();
        if (GravityEngine.Instance().units == GravityScaler.Units.DIMENSIONLESS) {
            phyPosition = transform.position;
        } else {
            phyPosition = nbody.initialPhysPosition;
        }
    }

    public bool IsFixed() {
		return true;
	}

	public void PreEvolve(float physicalScale, float massScale) {
		// nothing to do
	}
	
	public void Evolve(double physicsTime, float physicalScale, ref float[] r) {
        // dynamic origin shifting may move the position around
        r[0] = phyPosition.x;
        r[1] = phyPosition.y;
        r[2] = phyPosition.z;
	}

    public Vector3 GetVelocity() {
        return Vector3.zero;
    }

    public Vector3 GetPosition() {
        return phyPosition;
    }

    public void Move(Vector3 moveBy) {
        phyPosition += moveBy;
    }

    public void GEUpdate(GravityEngine ge) {
        // MapToScene may change things,so need to map every frame
        transform.position = ge.MapToScene(phyPosition);
    }
}
