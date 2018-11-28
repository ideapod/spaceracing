using System.Collections;
using System.Collections.Generic;
using UnityEngine;

/// <summary>
/// Launch to Orbit Mini-Game control logic
/// 
/// Provides user interface and game control for the launch and staging of a two stage rocket.
/// 
/// WASD orientation control for the ship
/// 
/// SPACE to launch/stage
/// 
/// F1/F2 for camera selection
/// 
/// Some UI is delegated:
/// - TimeZoom using 1-5 handled by the TmeZoom script
///-  Camera controls handled by CameraSpin script
///
/// </summary>
[RequireComponent(typeof(LineScaler))]
public class LaunchUI : MonoBehaviour
{

    public GameObject mainCameraBoom;
    public GameObject shipCameraBoom;

    //! Ship GO includes stages at launch
    public GameObject ship;

    // individual stages (are made independent objects upon staging)
    public GameObject[] stages;

    //! attached to spaceship model (model is a child of spceship NBody)
    public SpaceshipRV shipController;
   
    public float spinRate;

    // orbit parameter HUD - enabled at altitude HUD_ENABLE_ALTITUTE
    public GameObject orbitHUD;
    private const float HUD_ENABLE_ALTITUTE = 50f; // km

    public ShipInfo shipInfo;

    private bool running;
    private NBody shipNbody;
    private MultiStageEngine multiStage;

    private float launchTime;

    private const float LAUNCH_LOCKOUT = 5f;

    private LineScaler lineScaler;

    private const float SHIP_AT_EARTH_SCALE = 700f;
    private const float EARTH_LINE_SCALE = 300f;
    private const float SHIP_LINE_SCALE = 0.25f;

    private bool orbitPredictorOn = false;

    // awkward - flag to set line scale one frame after orbit predictor is turned on. Ick.
    private bool doLineScale = false; 
 
    // Use this for initialization
    void Start() {
        mainCameraBoom.SetActive(false);
        shipCameraBoom.SetActive(true);
        shipNbody = ship.GetComponent<NBody>();
        multiStage = ship.GetComponent<MultiStageEngine>();
        lineScaler = GetComponent<LineScaler>();

    }


    private void SelectCamera() {
        if (Input.GetKeyUp(KeyCode.F1)) {
            // Earth cam
            mainCameraBoom.SetActive(true);
            shipCameraBoom.SetActive(false);
            lineScaler.SetZoom(EARTH_LINE_SCALE);
            shipController.transform.localScale = SHIP_AT_EARTH_SCALE * Vector3.one;
            // Need to move forward so not hidden under the earth?
            
        } else if (Input.GetKeyUp(KeyCode.F2)) {
            // ship cam
            mainCameraBoom.SetActive(false);
            shipCameraBoom.SetActive(true);
            lineScaler.SetZoom(SHIP_LINE_SCALE);
            shipController.transform.localScale = Vector3.one;
        }
    }

    // J to stop engine and stage
    private void EngineControls() {
        if (Input.GetKeyUp(KeyCode.J)) {
            multiStage.NextStage();
            GravityEngine.Instance().TrajectoryRestart();
        }
        if (Input.GetKeyUp(KeyCode.X)) {
            multiStage.SetEngine(false);
        }
    }

    private bool UpdateThrustAxis() {
        // check for rotations
        Quaternion rotation = Quaternion.identity;
        // ship controller will rotate the engine thrust
        bool updated = false;
        if (Input.GetKeyUp(KeyCode.A)) {
            shipController.Rotate(spinRate, Vector3.forward);
            updated = true;
        } else if (Input.GetKeyUp(KeyCode.D)) {
            shipController.Rotate(-spinRate, Vector3.forward);
            updated = true;
        } else if (Input.GetKeyUp(KeyCode.W)) {
            shipController.Rotate(spinRate, Vector3.right);
            updated = true;
        } else if (Input.GetKeyUp(KeyCode.S)) {
            shipController.Rotate(-spinRate, Vector3.right);
        }
        return updated;
    }

    /// <summary>
    /// Detach the active stage and enable the next one if present. 
    /// </summary>
    private void DropStage() {
        // automatic staging when a preceeding stage runs out of fuel (engine still on)
        if (multiStage.activeStage < multiStage.numStages) {
            // release stage as child, add to GravityEngine
            GameObject spentStage = stages[multiStage.activeStage];
            spentStage.transform.parent = null;
            spentStage.AddComponent<NBody>();
            NBody nbody = spentStage.GetComponent<NBody>();
            nbody.initialPos = spentStage.transform.position;
            nbody.vel = shipNbody.vel;
            GravityEngine.Instance().AddBody(spentStage);

            // activate the next stage
            multiStage.NextStage();
            multiStage.SetEngine(true);
            GravityEngine.Instance().TrajectoryRestart();
            DisplayManager.Instance().DisplayMessage("Staging");

        } else {
            multiStage.SetEngine(false);
            DisplayManager.Instance().DisplayMessage("Engine Shutdown");
        }

    }

    // Update is called once per frame
    void Update() {

        SelectCamera();
        if (Input.GetKeyUp(KeyCode.Space)) {
            if (!running) {
                // Launch:
                // add the ship to the NBody integrator. The ship has a RocketEngine that 
                // will start immediatly. 
                // GE is already running so that e.g. objects in orbit are moving
                // Test - allow GE to start on launch
                if (!GravityEngine.Instance().evolveAtStart) {
                    GravityEngine.Instance().SetEvolve(true);
                }

                GravityEngine.Instance().AddBody(ship);
                running = true;
                multiStage.SetEngine(true);
                launchTime = GravityEngine.Instance().GetPhysicalTime();
                DisplayManager.Instance().DisplayMessage("Ignition");
            } else {
                DropStage();
            }
        }

        if (Input.GetKeyUp(KeyCode.P)) {
            orbitHUD.SetActive(true);
        }

        // toggle trajectory prediction
        if (Input.GetKeyUp(KeyCode.T)) {
            bool pred = GravityEngine.Instance().trajectoryPrediction; 
            GravityEngine.Instance().SetTrajectoryPrediction(!pred);
        }

         if (running) {
            EngineControls();
            bool updated = UpdateThrustAxis();
            float time = GravityEngine.Instance().GetPhysicalTime() - launchTime;
            float fuel = multiStage.GetFuel(GravityEngine.Instance().GetPhysicalTime());
            bool engineOn = multiStage.IsEngineOn();
            if (engineOn) {
                if (fuel > 0) {
                    if (updated) {
                        GravityEngine.Instance().TrajectoryRestart();
                    }
                } else {
                    DropStage();
                }
            }

            shipInfo.SetTextInfo(fuel, time, shipController.transform.rotation.eulerAngles);
            float altitude = shipInfo.GetAltitude();
            if ((altitude < 0) && (time > LAUNCH_LOCKOUT)) {
                // ship has crashed
                DisplayManager.Instance().DisplayMessage("Impact");
                // just remove the ship (so e.g. space station etc keep orbiting)
                GravityEngine.Instance().RemoveBody(ship);
            }
            // Awkward
            if(doLineScale) {
                lineScaler.FindAll();
                if (mainCameraBoom.activeInHierarchy) {
                    lineScaler.SetZoom(EARTH_LINE_SCALE);
                } else {
                    lineScaler.SetZoom(SHIP_LINE_SCALE);
                }

            }

            if ((altitude > HUD_ENABLE_ALTITUTE) && !orbitPredictorOn) {
                orbitHUD.SetActive(altitude > HUD_ENABLE_ALTITUTE);
                // have added an orbit predictor, need to apply correct scale, but LineScaler will not be there until next frame
                doLineScale = true; 
                orbitPredictorOn = true;
            }
        }
    }
}
