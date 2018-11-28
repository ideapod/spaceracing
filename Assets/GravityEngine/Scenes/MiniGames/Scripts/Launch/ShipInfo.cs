using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UI;


public class ShipInfo : MonoBehaviour {

    public Text time; 
    public Text position;
    public Text velocity;
    public Text altitude;
    public Text fuelText;
    public Text attitude;

    public NBody ship;

    private float altitudeNow;

    private const string time_format     =  "         Time: {0:0.0} ( x {1})";
    private const string pos_format      =  "     Pos [km]:   ({0:0.0}, {1:0.0}, {2:0.0})";
    private const string vel_format      =  "  Vel [km/hr]: ({0:0.0}, {1:0.0}, {2:0.0})";
    private const string fuel_format     =  "    Fuel (kg): {0:0.0}";
    private const string attitude_format =  "Pitch: {0:0.0} Roll: {1:0.0} Yaw: {2:0.0}";
    private const string altitude_format =  "Altitude (km): {0:0.0}";
    private float initial_altitude; 

    // Use this for initialization
    void Start () {
        GravityEngine ge = GravityEngine.Instance();
        time.text = string.Format(time_format,ge.GetPhysicalTime(), ge.GetTimeZoom() );
        Vector3 pos = ship.initialPos; 
        position.text = string.Format(pos_format, pos.x, pos.y, pos.z);
        velocity.text = string.Format(vel_format, 0f, 0f, 0f);
        RocketEngine rocketEngine = ship.GetComponent<RocketEngine>();
        if (rocketEngine == null)
        {
            Debug.LogError("Need a rocket component");
        }
        fuelText.text = string.Format(fuel_format, rocketEngine.GetFuel(ge.GetPhysicalTime()));
        attitude.text = string.Format(attitude_format, 0f, 0f, 0f);

        // assume earth at (0,0,0)
        initial_altitude = Vector3.Magnitude(ship.transform.position);
        altitude.text = string.Format(altitude_format, 0f);
    }

    public void SetTextInfo(float fuel, float timeValue, Vector3 angles)
    {
        time.text = string.Format(time_format, timeValue, GravityEngine.Instance().GetTimeZoom() );
        double[] r = new double[3];
        double[] v = new double[3];
        GravityEngine.Instance().GetPositionVelocityScaled(ship, ref r, ref v);

        position.text = string.Format(pos_format, r[0], r[1], r[2]);
        velocity.text = string.Format(vel_format, v[0], v[1], v[2]);
        fuelText.text = string.Format(fuel_format, fuel);
        attitude.text = string.Format(attitude_format, angles.y, angles.x, angles.z);

        altitudeNow = Vector3.Magnitude(ship.transform.position) - initial_altitude;
        altitude.text = string.Format(altitude_format, altitudeNow);
    }

    public float GetAltitude() {
        return altitudeNow;
    }


}
