using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class RotateAroundTarget : MonoBehaviour {

	public GameObject	target;
	public float 		distance;
	public float 		secondsPerRotation;
	float				rotationsPerSecond;

	// Use this for initialization
	void Start () {
		rotationsPerSecond = 1.0f / secondsPerRotation;
	}

	// Update is called once per frame
	void Update () {
		gameObject.transform.position = new Vector3(
			Mathf.Sin(2 * Mathf.PI * Time.fixedTime * rotationsPerSecond) * distance + target.transform.position.x,
			target.transform.position.y,
			Mathf.Cos(2 * Mathf.PI * Time.fixedTime * rotationsPerSecond) * distance + target.transform.position.z);

		transform.LookAt(target.transform.position);
	}
}
