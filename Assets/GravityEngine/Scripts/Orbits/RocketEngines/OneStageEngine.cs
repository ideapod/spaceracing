using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class OneStageEngine : RocketEngine {

    /// <summary>
    /// Implementation for a chemical rocket engine. Inputs are in the SI units
    /// indicated.  
    /// </summary>

    //! mass of empty ship (no fuel) in kg
    public double massEmpty;

    //! mass of fuel in kg
    public double massFuel;

    //! burn rate in kg/sec
    public double burnRate; 
    private double burnRateScaled; 

    //! thrust of engine in N
    public double thrust;

    private double tStart; 
    private double accelerationConversion;

    //! "live" thrust vector reflecintg orientation of ship model 
    private Vector3 thrustDirection;

    void Start()
    {
        SetEngine(engineOn);
        // Convert from real world units
        // Rocket engine thrust is given in SI units, so acceleration is in SI units
        // For Orbital scale need to convert to km/hr^2
        accelerationConversion = GravityScaler.AccelSItoGEUnits();
        // and convert to the physics scale used in the engine
        accelerationConversion = accelerationConversion * GravityScaler.AccelGEtoInternalUnits();

        burnRateScaled = burnRate * GravityScaler.GetGameSecondPerPhysicsSecond();
        thrustDirection = thrustAxis;
    }


    public override void SetEngine(bool on)
    {
        if (on)
        {
            tStart = GravityEngine.Instance().GetPhysicalTime();
        }
        engineOn = on;
    }

    public override double[] acceleration(double time)
    {
        // a(t) = Thrust/m(t) 
        // Will be called by both trajectory prediction and game evolution, so needs to be a function of time
        // (i.e. cannot reduce fuel each time this routine is called!)
        double[] a = new double[3] { 0, 0, 0 };
        if (engineOn)
        {
            double mass = massEmpty;
            double fuel = System.Math.Max(massFuel - burnRateScaled * (time - tStart), 0);
            if (fuel > 0)
            {
                double a_scalar = accelerationConversion * thrust / (mass + fuel);
                a[0] = a_scalar * thrustDirection.x;
                a[1] = a_scalar * thrustDirection.y;
                a[2] = a_scalar * thrustDirection.z;
            } 
            //Debug.Log(string.Format("a=({0}, {1}, {2}) fuel={3}", a[0], a[1], a[2], fuel));
        }

        return a;
    }

    public override float GetFuel(double time)
    {
        double fuel = System.Math.Max(massFuel - burnRateScaled * (time - tStart), 0);
        return (float) fuel;
    }

}
