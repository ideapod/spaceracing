
using UnityEngine;

/// <summary>
/// Rocket Engine Interface
/// Used to update acceleration in the massless body engine on a per-integration
/// step basis. 
/// 
/// Rocket engines are used for continuous application of force over a significant time
/// (unlike Maneuver which is used for an impulse change). 
/// 
/// Only one RocketEngine component can be added to an NBody and it MUST be added to a massless body
/// and have GE OptimizeMassless enabled. This means that if e.g. a multi-stage rocket with atmospheric drag
/// is required, the code must live in ONE RocketEngine implementation. (This is to eliminate traversing lists
/// etc. in the inner loop of the MasslessBodyEngine integrator)
/// 
/// It needs to provide the acceleration in the units used within the integrator:
/// - the choice of units in the GE is relevant
/// 
/// Concept:
/// By putting engine code deep in the integrator get several advantages:
/// - per integration step accuracy for rocket eqn
/// - ability to do trajectory prediction
/// - ability to timescale/timeZoom
/// 
/// The "cost" is that the engine implementation must not keep state - since the trajectory predictor will ask
/// it for it's mass at some future time. 
/// </summary>
public abstract class RocketEngine : MonoBehaviour {

    //! state of engine (or the active stage of a multi-stage engine)
    public bool engineOn;

    //! direction of thrust with respect to ship model. Must be unit length. 
    public Vector3 thrustAxis;

    //! "live" thrust vector reflecting orientation of ship model 
    protected Vector3 accelDirection;

    /// <summary>
    /// Determine the acceleration at the specified time. The accleration returned must be in the physical 
    /// units used by the integrators. For GE units other than DIMENSIONLESS this will involve some use of
    /// GravityScalar conversion functions. See the ChemicalRocket implementation for an example. 
    /// 
    /// When used in trajectory mode, this function may be called more than once for a given time, and
    /// as thrust changes. It will also be called for future times to determine the future path of the rocket.
    /// </summary>
    /// <param name="time">time for which the thrust should be determined.  </param>
    /// <returns></returns>
    public abstract double[] acceleration(double time);

    /// <summary>
    /// Set engine state to on. 
    /// </summary>
    /// <param name="state"></param>
    public abstract void SetEngine(bool state);

    public bool IsEngineOn() {
        return engineOn;
    }

    /// <summary>
    /// Get the fuel remaining at time T
    /// </summary>
    /// <param name="time"></param>
    /// <returns></returns>
    public abstract float GetFuel(double time);

    public void SetRotation(Quaternion rotation) {
        // TODO: Add gimbal angle changes
        accelDirection = rotation * (-1f * thrustAxis);
        Vector3 angles = rotation.eulerAngles;
        //Debug.LogFormat("Thrust direction:({0}, {1}, {2} ) angles=({3}, {4}, {5} ) thrustMag={6}", 
        //      accelDirection.x, accelDirection.y, accelDirection.z, angles.x, angles.y, angles.z, Vector3.Magnitude(accelDirection) );
        GravityEngine.Instance().TrajectoryRestart();
    }

}
