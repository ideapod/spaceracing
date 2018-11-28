using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class MultiStageEngine : RocketEngine
{

    /// <summary>
    /// Implementation for a multi-stage stage rocket with atmospheric drag that depends on velocity and altitude.
    /// 
    /// 
    /// </summary>


    public const int MAX_STAGES = 10;

    //! mass of payload (kg)
    public double payloadMass;

    // Per-stage details
    public int numStages;

    //! currently active stage - stages less than this are assumed jetisoned
    public int activeStage;

    //! Stage masses will be allocated by inspector script on first access
    // TODO - put all these parallel arrays into a struct-class

    //! mass of fuel in kg
    public double[] massFuel;

    //! mass of empty stage
    public double[] massStageEmpty;

    //! burn rate in kg/sec
    public double[] burnRate;

    //! burn rate scaled to physical units used in the integrator
    private double[] burnRateScaled;

    //! thrust of engine in N
    public double[] thrust;

    //! Particle effect object (flame etc.) for a stage. Set to active when engine is on. 
    public GameObject[] effectObject; 

#if UNITY_EDITOR
    // public for editor script control (persist foldout state etc.)
    public bool editorInited;
    public bool[] editorStageFoldout; 
#endif

    private double[] burnStart;
    private double accelerationConversion;


    void Start() {
        // Convert from real world units
        // Rocket engine thrust is given in SI units, so acceleration is in SI units
        // For Orbital scale need to convert to km/hr^2
        accelerationConversion = GravityScaler.AccelSItoGEUnits();
        // and convert to the physics scale used in the engine
        accelerationConversion = accelerationConversion * GravityScaler.AccelGEtoInternalUnits();

        burnStart = new double[MAX_STAGES];
        burnRateScaled = new double[MAX_STAGES];
        for (int i = 0; i < numStages; i++) {
            burnRateScaled[i] = burnRate[i] * GravityScaler.GetGameSecondPerPhysicsSecond();
        }
        // accelration will be opposite to thrust direction
        accelDirection = -thrustAxis;
        activeStage = 0;

        SetEngine(engineOn);
    }


    public override void SetEngine(bool on) {
        double time = GravityEngine.Instance().GetPhysicalTime();
        if (on) {
            burnStart[activeStage] = time;
            if (effectObject[activeStage] != null) {
                effectObject[activeStage].SetActive(true);
            }
        } else {
            if (effectObject[activeStage] != null) {
                effectObject[activeStage].SetActive(false);
            }
        }
        engineOn = on;
    }

 
 
    public override double[] acceleration(double time) {
        // a(t) = Thrust/m(t) 
        // Called on every integration timestep - so favour speed over beauty
        // Will be called by both trajectory prediction and game evolution, so needs to be a function of time
        // (i.e. cannot reduce fuel each time this routine is called!)
        double[] a = new double[3] { 0, 0, 0 };

        // always need the mass for drag calc (even if engine is off)
        double mass = payloadMass;
        // TODO - could optimize and precompute (mass per stage - fuel) in that stage
        // add stage masses
        for (int i = activeStage; i < numStages; i++) {
            mass += massStageEmpty[i];
        }
        // add fuel of unused stages
        for (int i = activeStage + 1; i < numStages; i++) {
            mass += massFuel[i];
        }
        // thrust
        if (engineOn) {
            double fuel = System.Math.Max(massFuel[activeStage] - burnRateScaled[activeStage] * (time - burnStart[activeStage]), 0);
            if (fuel > 0) {
                double a_scalar = accelerationConversion * thrust[activeStage] / (mass + fuel);
                a[0] = a_scalar * accelDirection.x;
                a[1] = a_scalar * accelDirection.y;
                a[2] = a_scalar * accelDirection.z;
            } 
        }

        return a;
    }

    /// <summary>
    /// Gets the fuel for the active stage
    /// </summary>
    /// <param name="time"></param>
    /// <returns></returns>
    public override float GetFuel(double time) {
        double fuel = massFuel[activeStage];
        if (engineOn) {
            fuel = System.Math.Max(massFuel[activeStage] - burnRateScaled[activeStage] * (time - burnStart[activeStage]), 0);
        }
        return (float)fuel;
    }

    /// <summary>
    /// Start the next stage
    /// </summary>
    /// <returns></returns>
    public void NextStage() {
        if (effectObject[activeStage] != null) {
            effectObject[activeStage].SetActive(false);
        }

        if (activeStage < numStages) {
            activeStage++;
            if (effectObject[activeStage] != null) {
                effectObject[activeStage].SetActive(true);
            }
        } else {
            engineOn = false;
        }
    }

}
