using System.Collections;
using System.Collections.Generic;
using UnityEngine;

// WORK IN PROGRESS

public class CircularPlaneChange : OrbitTransfer
{

    public CircularPlaneChange(OrbitData fromOrbit, OrbitData toOrbit) : base(fromOrbit, toOrbit)
    {

        name = "Circular Plane Change";

        // find  delta in inclination and Omega (chobotov p91)
        // TODO: quadrant issues
        // float planeChangeAngle = 0f;

        // angle to get to line of nodes
        // float nodesAngle = 0f;
        if (Mathf.Abs(fromOrbit.inclination) < TransferCalc.DELTA_INCL)
        {
            // planeChangeAngle = toOrbit.inclination;
        }
        else
        {
            // angles in Radians
            float i1 = Mathf.Deg2Rad * fromOrbit.inclination;
            float i2 = Mathf.Deg2Rad * toOrbit.inclination;
            float dOmega = Mathf.Deg2Rad * (toOrbit.omega_uc - toOrbit.omega_uc);
            float theta = Mathf.Acos(Mathf.Cos(i1) * Mathf.Cos(i2) +
                        Mathf.Sin(i1) * Mathf.Sin(i2) * Mathf.Cos(dOmega));

            float u1 = Mathf.Acos(Mathf.Cos(i1) * Mathf.Cos(theta) - Mathf.Cos(i2)) /
                            (Mathf.Sin(i1) * Mathf.Sin(theta));
            Debug.Log(string.Format("theta={0} degrees u1={1} degrees", theta*Mathf.Rad2Deg, u1*Mathf.Rad2Deg));

        }
    }

    private void Later(OrbitData fromOrbit, OrbitData toOrbit)
    { 
        // find velocity vector perpendicular to r for circular orbit
        double[] r_ship = new double[3];
        double[] v_ship = new double[3];
        GravityEngine ge = GravityEngine.Instance();

        ge.GetPositionVelocityScaled(fromOrbit.nbody, ref r_ship, ref v_ship);
        double[] r_center = new double[3];
        double[] v_center = new double[3];
        ge.GetPositionVelocityScaled(fromOrbit.centralMass, ref r_center, ref v_center);

        Vector3 pos_ship = new Vector3((float)r_ship[0], (float)r_ship[1], (float)r_ship[2]);
        Vector3 vel_ship = new Vector3((float)v_ship[0], (float)v_ship[1], (float)v_ship[2]);
        Vector3 pos_center = new Vector3((float)r_center[0], (float)r_center[1], (float)r_center[2]);
        Vector3 vel_center = new Vector3((float)v_center[0], (float)v_center[1], (float)v_center[2]);
        Vector3 r = pos_ship - pos_center;

        // want velocity relative to central mass (it could be moving)
        vel_ship = vel_ship - vel_center;

        // to get axis of orbit, can take r x v
        Vector3 axis = Vector3.Normalize(Vector3.Cross(pos_ship, vel_ship));
        // vis visa for circular orbit
        float mu = fromOrbit.centralMass.mass * ge.massScale;
        float v_mag = Mathf.Sqrt(mu / Vector3.Magnitude(r));
        // positive v is counter-clockwise
        Vector3 v_dir = Vector3.Normalize(Vector3.Cross(axis, r));
        Vector3 v_circular = v_mag * v_dir;

        Maneuver m1;
        m1 = new Maneuver();
        m1.nbody = fromOrbit.nbody;
        m1.mtype = Maneuver.Mtype.vector;
        m1.velChange = v_circular - vel_ship;
        m1.dV = Vector3.Magnitude(m1.velChange);
        Debug.Log(string.Format("v_ship={0} v_circular={1} axis={2}", vel_ship, v_circular, axis));
        m1.worldTime = GravityEngine.Instance().GetPhysicalTime();
        maneuvers.Add(m1);
    }

    public override string ToString()
    {
        return name;
    }
}
