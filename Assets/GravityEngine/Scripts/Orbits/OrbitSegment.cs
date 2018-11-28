﻿using UnityEngine;
using System.Collections;

/// <summary>
/// Orbit Segment.
/// An in-scene object that will determine the future orbit based on the current velocity and show
/// a segment of the orbit between the body and the to position. The short or long way between the
/// positions can be specified (in the case of an ellipse). 
/// 
/// The destination location can be provided by a reference to a gameObject or can be set by script
/// using the SetDestination() method. 
/// 
/// Depending on the velocity the orbit may be an ellipse or a hyperbola. This class
/// creates a delegate of each type to compute the orbital path. 
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
public class OrbitSegment: MonoBehaviour  {

    //! Number of points to be used in the line renderer for the orbit plot
	public int numPoints = 100; 
    //! The body for which an orbit is to be predicted
	public GameObject body;
    //! The object which body is in orbit around. 
	public GameObject centerBody;
    //! Script code will set the velocity explicitly. Do not retreive automatically
    public bool velocityFromScript = false;
    public bool shortPath = true;

    //! Location of the destination (optional, can use SetDestination() instead)
    public GameObject destination;
    // destination point in GE internal units
    private Vector3 destPoint;
    

    // velocity of body when set explicitly by script
    private Vector3 velocity;

	private NBody nbody;  
	private NBody aroundNBody; 
	private OrbitData orbitData; 
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

    /// <summary>
    /// Set the destination object (point on the orbit where the line rendering should stop). 
    /// 
    /// Set in GE internal physical units. 
    /// </summary>
    /// <param name=""></param>
    public void SetDestination(Vector3 d) {
        destPoint = d;
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

            // determine destination
            if (destination != null) {
                destPoint = GravityEngine.Instance().UnmapFromScene(destination.transform.position);
                destPoint = destination.transform.position/GravityEngine.Instance().physToWorldFactor;
            }

            if (velocityFromScript) {
                orbitData.SetOrbitForVelocity(nbody, aroundNBody, velocity);
            } else if (GravityEngine.Instance().GetEvolve()) {
                // OrbitData will update and retrieve velocities
                orbitData.SetOrbitForVelocity(nbody, aroundNBody);
            }

            // Now there is MapToRender MUST use physics and not transform position
            Vector3 centerPos = GravityEngine.Instance().GetPhysicsPosition(aroundNBody);
            Vector3 bodyPos = GravityEngine.Instance().GetPhysicsPosition(nbody);

            // Is the resulting orbit and ellipse or hyperbola?
            if (orbitData.ecc < 1f) {
                ellipseBase.InitFromOrbitData(orbitData);
                lineR.SetPositions(ellipseBase.OrbitSegmentPositions(numPoints, 
                                        centerPos, bodyPos, destPoint, shortPath));
            } else {
                orbitHyper.InitFromOrbitData(orbitData);
                bool mapToScene = true;
                lineR.SetPositions(orbitHyper.OrbitPositions(numPoints, centerPos, mapToScene));
            }
		}
	}
}