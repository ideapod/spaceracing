using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UI;

/// <summary>
/// Script to flip Camera parent to one of the anchors based on Function keys
/// </summary>
public class CameraAnchor : MonoBehaviour {

    public GameObject cameraBase; 
    public GameObject[] anchors;

    private Vector3 startPosition; 

	// Use this for initialization
	void Start () {
        startPosition = cameraBase.transform.position;
	}
	
	// Update is called once per frame
	void Update () {
		if (Input.GetKeyDown(KeyCode.F1)) {
            cameraBase.transform.parent = null;
            cameraBase.transform.position = startPosition; 
            return;
        }
        for (int i=0; i < anchors.Length; i++) {
            if (Input.GetKeyDown(KeyCode.F2+i)) {
                cameraBase.transform.parent = anchors[i].transform;
                cameraBase.transform.position = anchors[i].transform.position;
                return;
            }
        }
	}
}
