using System.Collections;
using System.Collections.Generic;
using UnityEngine;

/// <summary>
/// Gravity state.
/// Basic container to hold mass/position for massive bodies. 
///
/// </summary>
public class GravityState {

    public const int NDIM = 3; // Here to "de-magic" numbers. Some integrators have 3 baked in. Do not change.

    public int numBodies;

    //! masses of the massive bodies in the engine
    public double[] m;
    //! physics positions of massive objects in the engine
    public double[,] r;
    //! per GameObject flags for the integrator (INACTIVE, FIXED, TRAJ_DATA are current bits)
    public byte[] info;

    public List<GravityEngine.FixedBody> fixedBodies;

    // need to keep size^2 for simple collision detection between particles and 
    // massive bodies. Collisions between massive bodies are left to usual Unity
    // collider intrastructure
    public double[] size2; // size^2 used for detecting particle incursions

    //! size of the arrays (may exceed the number of bodies due to pre-allocation)
    public int arraySize;

    //! time of current state (in the engine physics time)
    public double time;

    //! State has trajectories that require updating
    public bool hasTrajectories;

    //!  physical time per evolver since start OR last timescale change
    public enum Evolvers { MASSIVE, MASSLESS, FIXED, PARTICLES };
    public double[] physicalTime;

    //! Flag to indicate running async (on a non-main thread). If true, cannot do debug logging or access scene
    public bool isAsync; 

    // Integrators - these are held in GE
    public INBodyIntegrator integrator;
    public MasslessBodyEngine masslessEngine;

    public List<GravityParticles> gravityParticles;

    // Delegate for handling Maneuvers. Access through GE wrapper methods, but
    // separate implementation in a delegate.
    public ManeuverMgr maneuverMgr;

    private const double EPSILON = 1E-4; 	// minimum distance for gravitatonal force

    private IForceDelegate forceDelegate;

    //! A force may be selective. Selective force needs to track integrator internal index structure
    private SelectiveForceBase selectiveForce; 

    private bool isCopy; 

    /// <summary>
    /// New Gravity state. Also need to call SetAlgorithmAndForce to fully configure. 
    /// </summary>
    /// <param name="size"></param>
    public GravityState(int size) {

        InitArrays(size);
        gravityParticles = new List<GravityParticles>();
        fixedBodies = new List<GravityEngine.FixedBody>();
        maneuverMgr = new ManeuverMgr();

#pragma warning disable 162     // disable unreachable code warning
        if (GravityEngine.DEBUG)
            Debug.Log("Created new (empty) gravityState");
#pragma warning restore 162

    }

    /// <summary>
    /// Clone constructor
    /// 
    /// Creates a deep copy suitable for independent evolution as a trajectory or for maneuver iterations. 
    /// Maneuvers will be executed but the callback to motify the owner of the maneuver will be skipped (only
    /// the real evolution will notify).
    /// </summary>
    /// <param name="fromState"></param>
    public GravityState(GravityState fromState) {
        m = new double[fromState.arraySize];
        r = new double[fromState.arraySize, NDIM];
        info = new byte[fromState.arraySize];
        size2 = new double[fromState.arraySize];
        physicalTime = new double[System.Enum.GetNames(typeof(Evolvers)).Length];
        arraySize = fromState.arraySize;
        numBodies = fromState.numBodies;

        // omitting hasTrajectories
        integrator = fromState.integrator.DeepClone();

        // don't copy particles, but need to init list
        gravityParticles = new List<GravityParticles>();

        // DO copy the maneuvers
        maneuverMgr = new ManeuverMgr(fromState.maneuverMgr);

        fixedBodies = new List<GravityEngine.FixedBody>(fromState.fixedBodies);

        for (int i = 0; i < physicalTime.Length; i++) {
            physicalTime[i] = fromState.physicalTime[i];
        }
        time = fromState.time;
        forceDelegate = fromState.forceDelegate;
        selectiveForce = fromState.selectiveForce;

        if (fromState.masslessEngine != null) {
            masslessEngine = fromState.masslessEngine.DeepClone();
            masslessEngine.ResetTrajectories((float)fromState.time);
        }

        for (int i = 0; i < arraySize; i++) {
            m[i] = fromState.m[i];
            info[i] = fromState.info[i];
            size2[i] = fromState.size2[i];
            r[i, 0] = fromState.r[i, 0];
            r[i, 1] = fromState.r[i, 1];
            r[i, 2] = fromState.r[i, 2];
        }

        // copies do not notify maneuver owners of maneuver completion. They are assumed to be "what if"
        // evolutions
        isCopy = true; 

#pragma warning disable 162     // disable unreachable code warning
        if (GravityEngine.DEBUG)
            Debug.Log("Created new (copy) gravityState");
#pragma warning restore 162
    }

    /// <summary>
    /// Set the integrator required for the chosen algorithm
    /// </summary>
    /// <param name="algorithm"></param>
    public void SetAlgorithmAndForce(GravityEngine.Algorithm algorithm, IForceDelegate forceDelegate) {
        this.forceDelegate = forceDelegate;
        // cast may be null if no force selection
        if (forceDelegate is SelectiveForceBase) {
            selectiveForce = (SelectiveForceBase)forceDelegate;
        }
        switch (algorithm) {
            case GravityEngine.Algorithm.LEAPFROG:
                integrator = new LeapFrogIntegrator(forceDelegate);
                break;
            //case GravityEngine.Algorithm.LEAPFROG_JOB:
            //    integrator = new LeapFrogJob(forceDelegate);
            //    break;
            case GravityEngine.Algorithm.HERMITE8:
                integrator = new HermiteIntegrator(forceDelegate);
                break;
            case GravityEngine.Algorithm.AZTRIPLE:
                integrator = new AZTripleIntegrator();
                break;
            default:
                Debug.LogError("Unknown algortithm");
                break;
        }
    }

    public void InitArrays(int arraySize) {
        m = new double[arraySize];
        r = new double[arraySize, NDIM];
        info = new byte[arraySize];
        size2 = new double[arraySize];

        physicalTime = new double[System.Enum.GetNames(typeof(Evolvers)).Length];
        this.arraySize = arraySize;
        if (selectiveForce) {
            selectiveForce.Init(arraySize);
        }
    }

    public bool GrowArrays(int growBy) {

        integrator.GrowArrays(growBy);

        double[] m_copy = new double[arraySize];
        double[,] r_copy = new double[arraySize, NDIM];
        byte[] info_copy = new byte[arraySize];
        double[] size2_copy = new double[arraySize];

        for (int i = 0; i < arraySize; i++) {
            m_copy[i] = m[i];
            r_copy[i, 0] = r[i, 0];
            r_copy[i, 1] = r[i, 1];
            r_copy[i, 2] = r[i, 2];
            info_copy[i] = info[i];
            size2_copy[i] = size2[i];
        }

        int newSize = arraySize + growBy;
        m = new double[newSize];
        r = new double[newSize, NDIM];
        info = new byte[newSize];
        size2 = new double[newSize];

        for (int i = 0; i < arraySize; i++) {
            m[i] = m_copy[i];
            info[i] = info_copy[i];
            size2[i] = size2_copy[i];
            r[i, 0] = r_copy[i, 0];
            r[i, 1] = r_copy[i, 1];
            r[i, 2] = r_copy[i, 2];
        }
        arraySize += growBy;

        if (selectiveForce) {
            selectiveForce.IncreaseToSize(arraySize);
        }

        return true;
    }

    public void Clear() {
        numBodies = 0;
        integrator.Clear();
        masslessEngine = null;
        gravityParticles.Clear();
        fixedBodies.Clear();
    }


    public void ResetPhysicalTime() {
        for (int i = 0; i < physicalTime.Length; i++) {
            physicalTime[i] = 0.0;
        }
    }

    public void AddFixedBody(GravityEngine.FixedBody fixedBody) {
        // Need to maintain order by kepler depth
        int insertAt = fixedBodies.Count;
        for (int i = 0; i < fixedBodies.Count; i++) {
            if (fixedBody.kepler_depth < fixedBodies[i].kepler_depth) {
                insertAt = i;
                break;
            }
        }
        fixedBodies.Insert(insertAt, fixedBody);
    }


    // Integrators need to use known positions to pre-determine accel. etc. to 
    // have valid starting values for evolution
    public void PreEvolve(GravityEngine ge) {
        // particles will pre-evolve when loading complete
        if (masslessEngine != null) {
            masslessEngine.PreEvolve(numBodies, this);
        }
        foreach (GravityEngine.FixedBody fixedBody in fixedBodies) {
            if (fixedBody.fixedOrbit != null) {
                fixedBody.fixedOrbit.PreEvolve(ge.physToWorldFactor, ge.massScale);
            }
        }
        integrator.PreEvolve(this, ref info);
    }

    /*******************************************
    * Main Physics Loop
    ********************************************/

    /// <summary>
    /// Evolve the objects subject to gravity. 
    /// 
    /// Normal evolution is done by passing in worldState with a time interval corresponding to the 
    /// frame advance time multiplied by the time zoom. 
    /// 
    /// For trajectory updates in the case where the trajectory is up to date, this will be for the same interval
    /// but starting at a future time. 
    /// 
    /// In order for trajectory computation to "catch up", there are times when the interval may be longer (but limited
    /// by the re-compute factor to avoid a huge recomputation on a single frame). 
    /// </summary>
    ///
    /// <param name="ge">The Gravity engine</param>
    /// <param name="physicsDt">The amount of physics DT to be evolved</param>
    /// 
    public bool Evolve(GravityEngine ge, double toPhysicsTime) {
        double gameDt = toPhysicsTime;
        bool trajectoryRestart = false;
        if (maneuverMgr.HaveManeuvers()) {
            List<Maneuver> maneuversInDt = maneuverMgr.ManeuversUntil((float) time);
            if (maneuversInDt.Count > 0) {
                foreach (Maneuver m in maneuversInDt) {
                    // evolve up to the time of the earliest maneuver
                    gameDt = System.Math.Max(m.worldTime -
                        (physicalTime[(int)GravityState.Evolvers.MASSIVE]), 0.0);
                    EvolveForTime(ge, gameDt);
                    m.Execute(this);
                    maneuverMgr.Executed(m, isCopy);
                }
                // recompute remaining time to evolve
                gameDt = System.Math.Max(time -
                        (physicalTime[(int)GravityState.Evolvers.MASSIVE]), 0.0);
                // if trajectories have made predictions, these need to be re-done since a manuever has
                // occured
                trajectoryRestart = true;  
            }
        }
        EvolveForTime(ge, gameDt);
        return trajectoryRestart;
    }

    private void EvolveForTime(GravityEngine ge, double physicsDt) {
        // Objective is to keep physical time proportional to game time 
        // Each integrator will run for at least as long as it is told but may overshoot
        // so correct time on next iteration. 
        // 
        // Keep the current physical time each integrator has reached in physicalTime[integrator_type]
        //
        double engineDt = ge.engineDt;
        if (physicsDt < engineDt)
            return;

        double timeEvolved = 0;

        // Need to move the integrators forward concurrently in steps matching the engineDt
        // - Hermite may be using a different timestep than this
        // - particles likely use a much longer timestep

        while (timeEvolved < physicsDt) {
            //==============================
            // Massive bodies
            //==============================
            // evolve all the massive game objects 
            double massiveDt = 0.0;
            massiveDt = integrator.Evolve(engineDt, this, ref info);
            physicalTime[(int) Evolvers.MASSIVE] += massiveDt;
            timeEvolved += massiveDt;
            // Debug.Log(string.Format("gameDt={0} integated={1} ptime={2} wtime={3}", gameDt, dt, physicalTime, worldTime));
            // LF is built in to particles and massless routines. They have their own DT built in
            // these run on a fixed timestep (if it is varied energy conservation is wrecked)
            // Track their evolution vs wall clock time seperately

            //==============================
            // Fixed Bodies
            //==============================
            // Update fixed update objects (if any)
            // Evolution is to a specific time - so use massive object physical time
            float[] r_new = new float[NDIM];
            foreach (GravityEngine.FixedBody fixedBody in fixedBodies) {
                fixedBody.fixedOrbit.Evolve(physicalTime[(int) Evolvers.MASSIVE],
                                            ge.physToWorldFactor,
                                            ref r_new);
                r[fixedBody.index, 0] = r_new[0];
                r[fixedBody.index, 1] = r_new[1];
                r[fixedBody.index, 2] = r_new[2];
            }

            //==============================
            // Particles (should only be present in worldState)
            //==============================
            if (gravityParticles.Count > 0) {
                double particle_dt = ge.GetParticleDt(); 
                if (physicalTime[(int)GravityState.Evolvers.PARTICLES] <
                        physicalTime[(int)GravityState.Evolvers.MASSIVE]) {
                    double evolvedFor = 0.0;
                    if (forceDelegate != null) {
                        foreach (GravityParticles nbp in gravityParticles) {
                            evolvedFor = nbp.EvolveWithForce(particle_dt, numBodies, this,
                                                        ref size2, ref info, forceDelegate);
                        }
                    } else {
                        foreach (GravityParticles nbp in gravityParticles) {
                            evolvedFor = nbp.Evolve(particle_dt, numBodies, this, ref size2, ref info);
                        }
                    }
                    physicalTime[(int) Evolvers.PARTICLES] += evolvedFor;
                }
            }

            //==============================
            // Massless
            //==============================
            if (masslessEngine != null) {
                // rockets need the time
                if (physicalTime[(int)  Evolvers.MASSLESS] <
                        physicalTime[(int) Evolvers.MASSIVE]) {
                    if (forceDelegate != null) {
                        physicalTime[(int) Evolvers.MASSLESS] +=
                                masslessEngine.EvolveWithForce(engineDt, time, numBodies, this, ref info, forceDelegate);
                    } else {
                        physicalTime[(int)Evolvers.MASSLESS] +=
                                masslessEngine.Evolve(engineDt, time, numBodies, this, ref info);
                    }
                }
            }
            // must update time so trajectory times are up to date
            time = physicalTime[(int)Evolvers.MASSIVE];
            if (hasTrajectories) {
                ge.UpdateTrajectories();
            }

        } // while
    }

    /// <summary>
    /// Remove the body at index and shuffle up the rest. Ensure the integrator does the same to stay
    /// in alignment. 
    /// </summary>
    /// <param name="index"></param>
    public void RemoveBodyAtIndex(int index) {
        integrator.RemoveBodyAtIndex(index);
        // shuffle the rest down, update indices
        for (int j = index; j < (numBodies - 1); j++) {
            info[j] =info[j + 1];
            m[j] = m[j + 1];
            r[j, 0] =r[j + 1, 0];
            r[j, 1] =r[j + 1, 1];
            r[j, 2] =r[j + 1, 2];
        }
        numBodies--;
        if (selectiveForce) {
            selectiveForce.RemoveBody(index);
        }
    }

    /// <summary>
    /// Get the internal position used by the physics engine. 
    /// </summary>
    /// <param name="nbody"></param>
    /// <returns></returns>
    public Vector3 GetPhysicsPosition(NBody nbody) {

        if (nbody == null || nbody.engineRef == null) {
            // may occur due to startup sequencing
            return Vector3.zero;
        }
        if (nbody.engineRef.bodyType == GravityEngine.BodyType.MASSLESS) {
            return masslessEngine.GetPosition(nbody);
        }
        if ((nbody.engineRef.fixedBody != null) && (nbody.engineRef.fixedBody.fixedOrbit != null)) {
            // If Kepler evolution 
            return nbody.engineRef.fixedBody.fixedOrbit.GetPosition();
        }
        return new Vector3((float) r[nbody.engineRef.index, 0],
                           (float) r[nbody.engineRef.index, 1],
                           (float) r[nbody.engineRef.index, 2]);
    }

    /// <summary>
    /// Get the physics velocity for an NBody as a double[]. 
    /// 
    /// </summary>
    /// <param name="nbody"></param>
    /// <param name="vel"></param>
    public void GetVelocityDouble(NBody nbody, ref double[] vel) {
        if ((nbody.engineRef.fixedBody != null) && (nbody.engineRef.fixedBody.fixedOrbit != null)) {
            // If Kepler evolution 
            Vector3 v = nbody.engineRef.fixedBody.fixedOrbit.GetVelocity();
            vel[0] = v.x;
            vel[1] = v.y;
            vel[2] = v.z;

        } else if (nbody.engineRef.bodyType == GravityEngine.BodyType.MASSLESS) {
            masslessEngine.GetVelocityDouble(nbody, ref vel);
        } else {
            integrator.GetVelocityDoubleForIndex(nbody.engineRef.index, ref vel);
        }
    }

    /// <summary>
    /// Set the physics velocity from a double array. 
    /// </summary>
    /// <param name="nbody"></param>
    /// <param name="velocity"></param>
    public void SetVelocityDouble(NBody nbody, ref double[] velocity) {
        if (nbody.engineRef.bodyType == GravityEngine.BodyType.MASSLESS) {
            masslessEngine.SetVelocityAtIndexDouble(nbody.engineRef.index, ref velocity);
        } else {
            integrator.SetVelocityDoubleForIndex(nbody.engineRef.index, ref velocity);
        }
    }

    /// <summary>
    /// Return the internal physics engine mass. 
    /// </summary>
    /// <param name="nbody"></param>
    /// <returns></returns>
    public double GetMass(NBody nbody) {
        if (nbody.engineRef.bodyType == GravityEngine.BodyType.MASSLESS) {
            return 0;
        }
        return m[nbody.engineRef.index];
    }

}
