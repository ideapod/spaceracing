using System.Collections;
using System.Collections.Generic;
using UnityEngine;

/// <summary>
/// Calculate the universal Lambert transfer. This algorithm allows the transfer to be elliptic or hyperbolic 
/// as necessary to meet the desinated transfer time. 
/// 
/// This algorithm does not provide an answer when the transfer is 180 degrees 
/// (Battin is a good alternative in this case). 
/// 
/// This algorithm is more general than the Lambert minimum energy solution since it provides the path for
/// a specified time of flight and also the final velocity.
/// 
/// Code adapted from Vallado, Fundamentals of Astrophysics and Applications
/// https://celestrak.com/software/vallado-sw.asp
/// </summary>
public class LambertUniversal: OrbitTransfer
{
    public double[] v1;
    public double[] v2;

    private Vector3d r1;
    private Vector3d r2;

    private OrbitData innerOrbit;
    private OrbitData outerOrbit;

    //! flag indicating if ship (fromOrbit) is transfering to an outer orbit
    private bool innerToOuter; 

    private double mu;

    private double a_min;
    // private double e_min;
    private double t_min;

    // center position of the orbit. r1 and r2 are with respect to this position
    private Vector3d center3d;

    /// <summary>
    /// Setup for Lambert Universal construction. Calculation is done via ComputeXfer and may be 
    /// done more than once for different transit times. 
    /// 
    /// Initial conditions are passed in when created. 
    /// 
    /// The transfer assumes there is an active NBody at the fromOrbit to be maneuvered to the 
    /// orbit specified by the toOrbit. One of the from/to orbit data elements must have a 
    /// centerNBody. 
    /// 
    /// </summary>
    /// <param name="fromOrbit"></param>
    /// <param name="toOrbit"></param>
    /// <param name="shortPath">Take shortest path along the ellipse (if xfer is an ellipse)</param>

    public LambertUniversal(OrbitData _fromOrbit, OrbitData _toOrbit, bool shortPath) : base(_fromOrbit, _toOrbit) {

        name = "LambertUniversal";

        // Fundamentals of Astrodynamics and Applications, Vallado, 4th Ed., Algorithm 56 p475 
        // Take r0, r => a_min, e_min, t_min, v0

        double[] r1_array = new double[] { 0, 0, 0 };
        double[] r2_array = new double[] { 0, 0, 0 };
        double[] center = new double[] { 0, 0, 0 };

        GravityEngine ge = GravityEngine.Instance();

        // If Nbody is available get position directly. If not (target marker) then 
        // use orbit phase to compute position
        ge.GetPositionDouble(centerBody, ref center);
        center3d = new Vector3d(ref center);

        if (fromOrbit.nbody != null) {
            ge.GetPositionDouble(fromOrbit.nbody, ref r1_array);
        } else {
            Vector3 pos = toOrbit.GetPhysicsPositionforEllipse(fromOrbit.phase);
            r1_array = new double[] { pos.x, pos.y, pos.z };
        }

        if (toOrbit.nbody != null) {
            ge.GetPositionDouble(toOrbit.nbody, ref r2_array);
        } else {
            Vector3 pos = toOrbit.GetPhysicsPositionforEllipse(toOrbit.phase);
            r2_array = new double[] { pos.x, pos.y, pos.z };
        }

        r1 = new Vector3d(ref r1_array) - center3d;
        r2 = new Vector3d(ref r2_array) - center3d;
        mu = fromOrbit.centralMass.mass;

        if (fromOrbit.a < toOrbit.a) {
            innerOrbit = fromOrbit;
            outerOrbit = toOrbit;
            innerToOuter = true;
        } else {
            innerOrbit = toOrbit;
            outerOrbit = fromOrbit;
            innerToOuter = false;
        }
        mu = ge.GetPhysicsMass(fromOrbit.centralMass);
        // determine t_min. 
        // If Nbody is available get position directly. If not (target marker) then 
        // use orbit phase to compute position
        double[] r0 = new double[] { 0, 0, 0 };
        double[] r = new double[] { 0, 0, 0 };
        if (innerOrbit.nbody != null) {
            ge.GetPositionDouble(innerOrbit.nbody, ref r0);
        } else {
            Vector3 pos = toOrbit.GetPhysicsPositionforEllipse(innerOrbit.phase);
            r0 = new double[] { pos.x, pos.y, pos.z };
        }

        if (outerOrbit.nbody != null) {
            ge.GetPositionDouble(outerOrbit.nbody, ref r);
        } else {
            Vector3 pos = toOrbit.GetPhysicsPositionforEllipse(outerOrbit.phase);
            r = new double[] { pos.x, pos.y, pos.z };
        }

        ComputeMinTime(ref r0, ref r, fromOrbit.centralMass, shortPath);
    }

    /// <summary>
    /// Create a Lambert transfer from a given NBody orbit to a position in physics space. 
    /// 
    /// The constructor will compute the minimum energy transfer to the point taking the shortest
    /// path setting into account. The time for the transfer can then be retreived. 
    /// </summary>
    /// <param name="_fromOrbit"></param>
    /// <param name="r_from"></param>
    /// <param name="r_to"></param>
    /// <param name="shortPath"></param>
    public LambertUniversal(OrbitData _fromOrbit, Vector3d r_from, Vector3d r_to, bool shortPath) : base(_fromOrbit) {
        name = "LambertUniversal";

        GravityEngine ge = GravityEngine.Instance();

        // If Nbody is available get position directly. If not (target marker) then 
        // use orbit phase to compute position
        center3d = ge.GetPositionDoubleV3(centerBody);
        // TODO: require r1 < r2
        r1 = r_from - center3d;
        r2 = r_to - center3d;

        fromOrbit = _fromOrbit;
        innerOrbit = fromOrbit;
        innerToOuter = true;
        mu = fromOrbit.centralMass.mass;

        mu = ge.GetPhysicsMass(fromOrbit.centralMass);
        // determine t_min. 
        double[] r0 = new double[] { r1.x, r1.y, r1.z };
        double[] r = new double[] { r2.x, r2.y, r2.z };
        ComputeMinTime(ref r0, ref r, _fromOrbit.centralMass, shortPath);

    }

    /// <summary>
    /// Compute the transfer time for the minimum energy transfer. 
    /// </summary>
    /// <returns></returns>
    private double ComputeMinTime(ref double[] r0, ref double[] r, NBody centerBody, bool shortPath) {

        // Fundamentals of Astrodynamics and Applications, Vallado, 4th Ed., Algorithm 56 p475 
        // Take r0, r => a_min, e_min, t_min
        // This approach assumes r > r0. 
        // If we have r0 > r then need to flip them and then reverse the velocities when doing the 
        // maneuver determination

        double[] center = new double[] { 0, 0, 0 };

        GravityEngine ge = GravityEngine.Instance();
        ge.GetPositionDouble(centerBody, ref center);
        Vector3d center3d = new Vector3d(ref center);

        Vector3d r0_vec = new Vector3d(ref r0) - center3d;
        Vector3d r_vec = new Vector3d(ref r) - center3d;

        double r0_mag = Vector3d.Magnitude(r0_vec);
        double r_mag = Vector3d.Magnitude(r_vec);
        double cos_delta_nu = Vector3d.Dot(r0_vec, r_vec) / (r0_mag * r_mag);

        double c = System.Math.Sqrt(r0_mag * r0_mag + r_mag * r_mag
                    - 2 * r0_mag * r_mag * cos_delta_nu);
        double s = (r0_mag + r_mag + c) / 2;
        a_min = s / 2;
        // double p_min = r0_mag * r_mag * (1 - cos_delta_nu) / c;
        // e_min = System.Math.Sqrt(1 - 2 * p_min / s);

        // alpha for elliptical xfer, min energy
        double alpha_e = System.Math.PI;

        // p471: set beta_e negative if delta_nu > 18
        double beta_e = 2.0 * System.Math.Asin(System.Math.Sqrt((s - c) / s));
        // t_min for a_min
        // +/- in front of (beta_e...). Which sign? Negative for shorter path. 
        double sign = -1;
        if (!shortPath) {
            sign = 1;
        }
        t_min = System.Math.Sqrt(a_min * a_min * a_min / mu) * (alpha_e + sign * (beta_e - System.Math.Sin(beta_e)));
        return t_min;
    }

    // only need value up to 13 in the Lambert universal, so make it a look up
    private double[] factorial = { 1.0, 1.0, 2.0, 6.0, 24.0, 120.0, 720.0, 5040.0, 40320.0,
            362880.0, 3628800.0, 3.991680E7, 4.790016E8, 6227021E9};

    private double Factorial(int x) {
        return factorial[x];
    }

    private void FindC2C3(double znew, out double c2new, out double c3new) {
        double small, sqrtz;
        small = 0.00000001;

        // -------------------------  implementation   -----------------
        if (znew > small) {
            sqrtz = System.Math.Sqrt(znew);
            c2new = (1.0 - System.Math.Cos(sqrtz)) / znew;
            c3new = (sqrtz - System.Math.Sin(sqrtz)) / (sqrtz * sqrtz * sqrtz);
        } else {
            if (znew < -small) {
                sqrtz = System.Math.Sqrt(-znew);
                c2new = (1.0 - System.Math.Cosh(sqrtz)) / znew;
                c3new = (System.Math.Sinh(sqrtz) - sqrtz) / (sqrtz * sqrtz * sqrtz);
            } else {
                c2new = 0.5;
                c3new = 1.0 / 6.0;
            }
        }
    }  // findc2c3

    /// <summary>
    /// 
    /// </summary>
    /// <param name="reverse">direction of motion, true for reverse path (long way around)</param>
    /// <param name="df">Controls intial guess, undocumented in Vallado</param>
    /// <param name="nrev">Number of revolutions until transfer (typically 0)</param>
    /// <param name="dtsec">Required time of flight</param>
    /// <returns>error ()</returns>
    ///     error = 1;   // g not converged
    ///     error = 2;   // y negative
    ///     error = 3;   // impossible 180 transfer
    public int ComputeXfer(
        bool reverse,
        bool df,
        int nrev,
        double dtsec) 
    { 
        const double small = 0.0000001;
        const int numiter = 40;
        // const double mu = 398600.4418;  // m3s2
        const double pi = System.Math.PI;

        v1 = new double[] { 0, 0, 0 };
        v2 = new double[] { 0, 0, 0 };

        int loops, ynegktr;
        double vara, y, upper, lower, cosdeltanu, f, g, gdot, xold, xoldcubed, magr1, magr2,
            psiold, psinew, c2new, c3new, dtnew, c2dot, c3dot, dtdpsi, psiold2;

        y = 0.0; 

        /* --------------------  initialize values   -------------------- */
        int error = 0;
        psinew = 0.0;

        magr1 = r1.magnitude;
        magr2 = r2.magnitude;

        cosdeltanu = Vector3d.Dot(r1, r2) / (magr1 * magr2);
        if (reverse)
            vara = -System.Math.Sqrt(magr1 * magr2 * (1.0 + cosdeltanu));
        else
            vara = System.Math.Sqrt(magr1 * magr2 * (1.0 + cosdeltanu));

        /* -------- set up initial bounds for the bissection ------------ */
        if (nrev == 0) {
            upper = 4.0 * pi * pi;  // could be negative infinity for all cases
            lower = -4.0 * pi * pi; // allow hyperbolic and parabolic solutions
        } else {
            lower = 4.0 * nrev * nrev * pi * pi;
            upper = 4.0 * (nrev + 1.0) * (nrev + 1.0) * pi * pi;
        }

        /* ----------------  form initial guesses   --------------------- */
        psinew = 0.0;
        xold = 0.0;
        if (nrev == 0) {
            // use log to get initial guess
            // empirical relation here from 10000 random draws
            // 10000 cases up to 85000 dtsec  0.11604050x + 9.69546575
            psiold = (System.Math.Log(dtsec) - 9.61202327) / 0.10918231;
            if (psiold > upper)
                psiold = upper - pi;
        } else {
            if (df)
                psiold = lower + (upper - lower) * 0.3;
            else
                psiold = lower + (upper - lower) * 0.6;
        }
        FindC2C3(psiold, out c2new, out c3new);


        /* --------  determine if the orbit is possible at all ---------- */
        if (System.Math.Abs(vara) > small)  // 0.2??
        {
            loops = 0;
            ynegktr = 1; // y neg ktr
            dtnew = -10.0;
            while ((System.Math.Abs(dtnew - dtsec) >= small) && (loops < numiter) && (ynegktr <= 10)) {
                loops = loops + 1;
                if (System.Math.Abs(c2new) > small)
                    y = magr1 + magr2 - (vara * (1.0 - psiold * c3new) / System.Math.Sqrt(c2new));
                else
                    y = magr1 + magr2;
                /* ------- check for negative values of y ------- */
                if ((vara > 0.0) && (y < 0.0)) {
                    ynegktr = 1;
                    while ((y < 0.0) && (ynegktr < 10)) {
                        psinew = 0.8 * (1.0 / c3new) *
                            (1.0 - (magr1 + magr2) * System.Math.Sqrt(c2new) / vara);

                        /* ------ find c2 and c3 functions ------ */
                        FindC2C3(psinew, out c2new, out c3new);
                        psiold = psinew;
                        lower = psiold;
                        if (System.Math.Abs(c2new) > small)
                            y = magr1 + magr2 -
                            (vara * (1.0 - psiold * c3new) / System.Math.Sqrt(c2new));
                        else
                            y = magr1 + magr2;
                        ynegktr++;
                    }
                }

                if (ynegktr < 10) {
                    if (System.Math.Abs(c2new) > small)
                        xold = System.Math.Sqrt(y / c2new);
                    else
                        xold = 0.0;
                    xoldcubed = xold * xold * xold;
                    dtnew = (xoldcubed * c3new + vara * System.Math.Sqrt(y)) / System.Math.Sqrt(mu);

                    // try newton rhapson iteration to update psi
                    if (System.Math.Abs(psiold) > 1e-5) {
                        c2dot = 0.5 / psiold * (1.0 - psiold * c3new - 2.0 * c2new);
                        c3dot = 0.5 / psiold * (c2new - 3.0 * c3new);
                    } else {
                        psiold2 = psiold * psiold;
                        c2dot = -1.0 / Factorial(4) + 2.0 * psiold / Factorial(6) - 3.0 * psiold2 / Factorial(8) +
                            4.0 * psiold2 * psiold / Factorial(10) - 5.0 * psiold2 * psiold2 / Factorial(12);
                        c3dot = -1.0 / Factorial(5) + 2.0 * psiold / Factorial(7) - 3.0 * psiold2 / Factorial(9) +
                            4.0 * psiold2 * psiold / Factorial(11) - 5.0 * psiold2 * psiold2 / Factorial(13);
                    }
                    dtdpsi = (xoldcubed * (c3dot - 3.0 * c3new * c2dot / (2.0 * c2new)) + vara / 8.0 * (3.0 * c3new * System.Math.Sqrt(y) / c2new + vara / xold)) / System.Math.Sqrt(mu);
                    psinew = psiold - (dtnew - dtsec) / dtdpsi;

                    // check if newton guess for psi is outside bounds(too steep a slope)
                    if (System.Math.Abs(psinew) > upper || psinew < lower) {
                        // --------readjust upper and lower bounds------ -
                        if (dtnew < dtsec) {
                            if (psiold > lower)
                                lower = psiold;
                        }
                        if (dtnew > dtsec) {
                            if (psiold < upper)
                                upper = psiold;
                        }
                        psinew = (upper + lower) * 0.5;
                    }
                    /* -------------- find c2 and c3 functions ---------- */
                    FindC2C3(psinew, out c2new, out c3new);
                    psiold = psinew;

                    /* ---- make sure the first guess isn't too close --- */
                    if ((System.Math.Abs(dtnew - dtsec) < small) && (loops == 1))
                        dtnew = dtsec - 1.0;
                }
            }

            if ((loops >= numiter) || (ynegktr >= 10)) {
                error = 1; // g not converged

                if (ynegktr >= 10) {
                    error = 2;  // y negative
                }
            } else {
                /* ---- use f and g series to find velocity vectors ----- */
                f = 1.0 - y / magr1;
                gdot = 1.0 - y / magr2;
                g = 1.0 / (vara * System.Math.Sqrt(y / mu)); // 1 over g
                //	fdot = sqrt(y) * (-magr2 - magr1 + y) / (magr2 * magr1 * vara);
                //for (int i = 0; i < 3; i++) {
                //    v1[i] = ((r2[i] - f * r1[i]) * g);
                //    v2[i] = ((gdot * r2[i] - r1[i]) * g);
                //}
                v1[0] = ((r2.x - f * r1.x) * g);
                v2[0] = ((gdot * r2.x - r1.x) * g);
                v1[1] = ((r2.y - f * r1.y) * g);
                v2[1] = ((gdot * r2.y - r1.y) * g);
                v1[2] = ((r2.z - f * r1.z) * g);
                v2[2] = ((gdot * r2.z - r1.z) * g);
            }
        } else
            error = 3;   // impossible 180 transfer

        // Determine maneuvers needed for ship (fromOrbit)
        if (error == 0) {
            maneuvers.Clear();
            Vector3 v1_vec = new Vector3((float)v1[0], (float)v1[1], (float)v1[2]);
            Vector3 v2_vec = new Vector3((float)v2[0], (float)v2[1], (float)v2[2]);
            // Departure
            Maneuver departure = new Maneuver();
            departure.mtype = Maneuver.Mtype.vector;
            departure.nbody = fromOrbit.nbody;
            departure.physPosition = r1 + center3d;
            if (innerToOuter) { 
                departure.velChange = v1_vec - GravityEngine.Instance().GetVelocity(innerOrbit.nbody);
            } else {
                // Need to establish arrival velocity
                departure.velChange = v2_vec - GravityEngine.Instance().GetVelocity(outerOrbit.nbody);
            }
            departure.worldTime = (float)GravityEngine.Instance().GetWorldTime(); ;
            maneuvers.Add(departure);
            deltaV = departure.velChange.magnitude;

            // Arrival (will not be required if intercept)
            if (toOrbit != null) {
                Maneuver arrival = new Maneuver();
                arrival.nbody = fromOrbit.nbody;
                arrival.physPosition = r2 + center3d;
                arrival.worldTime = departure.worldTime + (float)dtsec;
                arrival.mtype = Maneuver.Mtype.vector;
                arrival.velChange = toOrbit.GetPhysicsVelocityForEllipse(toOrbit.phase) - v2_vec;
                maneuvers.Add(arrival);
                deltaV += arrival.velChange.magnitude;
            }
        }
        return error;
    }

    public Vector3 GetTransferVelocity() {
         return new Vector3((float)v1[0], (float)v1[1], (float)v1[2]);
    }

    public Vector3d GetTransferVelocityDouble() {
        return new Vector3d(v1[0],v1[1], v1[2]);
    }

    //! Time of flight for minimum energy trajectory
    public double GetTMin() {
        return t_min;
    }


}

