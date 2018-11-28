using UnityEngine;
using System.Collections;

/// <summary>
/// Orbit predictor.
/// An in-scene object that will determine the future orbit based on the current velocity. 
/// Depending on the velocity the orbit may be an ellipse or a hyperbola. This class
/// requires a delegate of each type to compute the orbital path. 
///
/// Orbit prediction is based on the two-body problem and is with respect to one other
/// body (presumably the dominant source of gravity for the affected object). The OrbitPredictor will
/// add both an OrbitEllipse and an OrbitHyper and use one or the other to plot the projected orbit
/// based on the velocity. The class OrbitData is used to determine the orbital parameters for the
/// velocity. 
///
/// The general orbit prediction problem is significantly harder - it requires simulating the
/// entire scene into the future - re-computing whenever user input is provided. This is provided by
/// the Trajectory prediction sub-system. 
/// </summary>
[RequireComponent(typeof(LineRenderer))]
public class OrbitPredictor: MonoBehaviour  {

    //! Number of points to be used in the line renderer for the orbit plot
	public int numPoints = 100; 
    //! The body for which an orbit is to be predicted
	public GameObject body;
    //! The object which body is in orbit around. 
	public GameObject centerBody;
    //! Script code will set the velocity explicitly. Do not retreive automatically
    public bool velocityFromScript = false;

    // velocity of body when set explicitly by script
    private Vector3 velocity;

	private NBody nbody;  
	private NBody aroundNBody; 
	public OrbitData orbitData; 
	private EllipseBase ellipseBase; 
	private OrbitHyper orbitHyper; 

	private LineRenderer lineR;

	// Use this for initialization
	void Start () {

		nbody = body.GetComponent<NBody>();
		if (nbody == null) {
			Debug.LogWarning("Cannot show orbit - Body requires NBody component");
			return;
		}
		lineR = GetComponent<LineRenderer>();
		lineR.positionCount = numPoints; 
		orbitData = new OrbitData(); 

		ellipseBase = transform.gameObject.AddComponent<EllipseBase>();
		orbitHyper = transform.gameObject.AddComponent<OrbitHyper>();

        if (centerBody == null) {
            Debug.LogError("No center body for OrbitPredictor");
        }
        SetCenterObject(centerBody);
	}

    public void SetCenterObject(GameObject newCenterBody) {
        if (centerBody != null) {
            centerBody = newCenterBody;
            aroundNBody = centerBody.GetComponent<NBody>();
            ellipseBase.centerObject = centerBody;
            orbitHyper.centerObject = centerBody;
        }
    }

    public void SetVelocity(Vector3 v) {
        velocity = v;
    }

    public Vector3 GetVelocity() {
        return velocity;
    }

    // Update is called once per frame
    void Update () {
		if ((nbody != null) && (centerBody != null)) {

            // velocities not normally updated in NBody (to reduce CPU). Need to ask
            // object to update its velocity from Gravity Engine if engine is running
            if (GravityEngine.Instance().GetEvolve()) {
                nbody.UpdateVelocity();
            }

            if (velocityFromScript) {
                orbitData.SetOrbitForVelocity(nbody, aroundNBody, velocity);
            } else if (GravityEngine.Instance().GetEvolve()) {
                // OrbitData will update and retrieve both velocities
                orbitData.SetOrbitForVelocity(nbody, aroundNBody);
            } else {
                // use Nbody velocity without updates
                orbitData.SetOrbitForVelocity(nbody, aroundNBody, nbody.vel_phys);
            }

            // Now there is MapToRender MUST use physics and not transform position
            Vector3 centerPos = GravityEngine.Instance().GetPhysicsPosition(aroundNBody);

            // Is the resulting orbit and ellipse or hyperbola?
            bool mapToScene = true; 
            if (orbitData.ecc < 1f) {
                ellipseBase.InitFromOrbitData(orbitData);
                lineR.SetPositions(ellipseBase.OrbitPositions(numPoints, centerPos, mapToScene));
            } else {
                orbitHyper.InitFromOrbitData(orbitData);
                lineR.SetPositions(orbitHyper.OrbitPositions(numPoints, centerPos, mapToScene));
            }
		}
	}
}
