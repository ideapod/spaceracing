using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class HohmannXfer : OrbitTransfer {

    private const float TWO_PI = 2f * Mathf.PI;

	public HohmannXfer(OrbitData fromOrbit, OrbitData toOrbit, bool rendezvous) : base(fromOrbit, toOrbit) {
		name = "Hohmann";
        if (rendezvous) {
            name += " Rendezvous";
        } else {
            name += " Transfer";
        }

		// Hohmann xfer is via an ellipse from one circle to another. The ellipse is uniquely
		// defined by the radius of from and to. 
		// Equations from Chobotov Ch 5.4
		float r_inner = 0f;
		float r_outer = 0f; 
		// From orbit result is in physics quantities
		if (fromOrbit.a < toOrbit.a) {
			r_inner = fromOrbit.a;
			r_outer = toOrbit.a;
		} else {
			r_inner = toOrbit.a;
			r_outer = fromOrbit.a;
		}

		// mass scale applied in Orbit data
		float v_inner = Mathf.Sqrt(fromOrbit.mu/r_inner);
		float rf_ri = r_outer/r_inner;
		float dV_inner = v_inner * (Mathf.Sqrt(2f*rf_ri/(1f+rf_ri)) -1f);

		float v_outer = Mathf.Sqrt(fromOrbit.mu/r_outer);
        //simplify per Roy (12.22)
		float dV_outer = v_outer * (1f - Mathf.Sqrt(2/(1+rf_ri)));

        // Debug.LogFormat("r_in={0} r_out={1}  v_inner={2} v_outer={3}", r_inner, r_outer, v_inner, v_outer);

		// transfer time
        // Need to flip rf_ri for inner orbits to get the correct transfer_time
        // (should re-derive for this case sometime to see why)
		float transfer_time = 0f;
        // time to wait for first maneuver (rendezvous case)
        float tWait = 0f;

        // Build the manuevers required
        deltaV = 0f;
		float worldTime = GravityEngine.Instance().GetPhysicalTime();

		Maneuver m1;
		m1 = new Maneuver();
		m1.nbody = fromOrbit.nbody;
		m1.mtype = Maneuver.Mtype.scalar;
		m1.worldTime = worldTime;
		Maneuver m2;
		m2 = new Maneuver();
		m2.nbody = fromOrbit.nbody;
		m2.mtype = Maneuver.Mtype.scalar;

        if (fromOrbit.a < toOrbit.a) {
            float subexpr = 1f + rf_ri;
            transfer_time = fromOrbit.period / Mathf.Sqrt(32) * Mathf.Sqrt(subexpr * subexpr * subexpr);
            if (rendezvous) {
                // need to determine wait time for first maneuver to phase the arrival
                // find angle by which outer body must lead (radians)
                // Chobotov 7.3
                float subexpr2 = 0.5f * (1f + r_inner / r_outer);
                float theta_h = Mathf.PI * (1f - Mathf.Sqrt(subexpr2 * subexpr2 * subexpr2));
                // find current angular seperation
                float phase_gap = toOrbit.phase - fromOrbit.phase;
                if (phase_gap < 0)
                    phase_gap += 360f;
                // need seperation to be theta_h 
                float dTheta = Mathf.Deg2Rad * phase_gap - theta_h;
                if (dTheta < 0) {
                    dTheta += TWO_PI;
                }
                // need to wait for phase_gap to reduce to this value. It reduces at a speed based on the difference 
                // in the angular velocities. 
                float dOmega = TWO_PI / fromOrbit.period - TWO_PI / toOrbit.period;
                tWait = dTheta / dOmega;
                //Debug.LogFormat("inner_phase= {0} out_phase={1} phase_gap(deg)={2} thetaH={3} dTheta(rad)={4} dOmega={5} tWait={6}", 
                //    fromOrbit.phase, toOrbit.phase, phase_gap, theta_h, dTheta, dOmega, tWait);

            }
            // from inner to outer
            // first maneuver is to a higher orbit
            m1.dV = dV_inner;
			deltaV += dV_inner;
			maneuvers.Add(m1);
			// second manuever is opposite to velocity
			m2.dV = dV_outer;
			deltaV += dV_outer;
			maneuvers.Add(m2);
		} else {
            float subexpr_in = 1f + r_inner/r_outer;
            transfer_time = fromOrbit.period / Mathf.Sqrt(32) * Mathf.Sqrt(subexpr_in * subexpr_in * subexpr_in);
            if (rendezvous) {
                // Chobotov 7.2/7.3 (modified for outer to inner, use (Pi+Theta) and Pf not Pi
                float subexpr2 = 0.5f * (1f + r_outer / r_inner);
                float theta_h = Mathf.PI * (1f + Mathf.Sqrt(subexpr2 * subexpr2 * subexpr2));
                // find current angular seperation
                float phase_gap = fromOrbit.phase - toOrbit.phase;
                if (phase_gap < 0)
                    phase_gap += 360f;
                // need seperation to be -theta_h 
                float dTheta = Mathf.Deg2Rad * phase_gap - theta_h;
                // Can need inner body to go around more than once...
                while (dTheta < 0) {
                    dTheta += TWO_PI;
                }
                // larger (inner) omega first 
                float dOmega =  TWO_PI / toOrbit.period - TWO_PI / fromOrbit.period;
                tWait = dTheta / dOmega;
                //Debug.LogFormat("inner_phase= {0} out_phase={1} phase_gap(deg)={2} thetaH(deg)={3} dTheta(rad)={4} dOmega={5} tWait={6}",
                //    toOrbit.phase, fromOrbit.phase, phase_gap, Mathf.Rad2Deg*theta_h, dTheta, dOmega, tWait);

            }
            // from outer to inner
            // first maneuver is to a lower orbit
            m1.dV = -dV_outer;
			deltaV += dV_outer;
			maneuvers.Add(m1);
			// second manuever is opposite to velocity
			m2.dV = -dV_inner;
			deltaV += dV_outer;
			maneuvers.Add(m2);
		}
        m1.worldTime = worldTime + tWait;
        m2.worldTime = worldTime + tWait + transfer_time;
        deltaT = tWait + transfer_time;

     }

    public override string ToString() {
		return name;
	}

    public HohmannXfer CreateTransferCopy(bool rendezvous) {

        HohmannXfer newXfer = new HohmannXfer(this.fromOrbit, this.toOrbit, rendezvous);
        return newXfer;
    }


}
