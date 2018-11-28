﻿using UnityEngine;
using UnityEditor;
using System.Collections;

[CustomEditor(typeof(OrbitPredictor), true)]
public class OrbitPredictorEditor : Editor {

	private const string centerTip = "NBody that the body is in orbit around.";
	private const string bodyTip = "NBody for orbit prediction.";
	private const string rTip = "Number of points to use in line renderering of orbit.";
    private const string vTip = "Velocity will be set explicitly by a script (do not ask GE for velocity every frame)";
    public override void OnInspectorGUI()
	{
		GUI.changed = false;
		OrbitPredictor orbit = (OrbitPredictor) target;

		GameObject body; 
		GameObject centerObject = orbit.centerBody;
        bool vFromScript = orbit.velocityFromScript;
		int numPoints;

        centerObject = (GameObject)EditorGUILayout.ObjectField(
                new GUIContent("CenterObject", centerTip),
                orbit.centerBody,
                typeof(GameObject),
                true);
 
        body = (GameObject) EditorGUILayout.ObjectField(
				new GUIContent("Body", bodyTip), 
				orbit.body,
				typeof(GameObject), 
				true);

        vFromScript = EditorGUILayout.Toggle(new GUIContent("Velocity From Script", vTip),  vFromScript);

		numPoints = EditorGUILayout.IntField(new GUIContent("Number of Points", rTip), orbit.numPoints);

		if (GUI.changed) {
			Undo.RecordObject(orbit, "OrbitEllipse Change");
			orbit.centerBody = centerObject;
			orbit.body = body;
			orbit.numPoints = numPoints;
            orbit.velocityFromScript = vFromScript;
			EditorUtility.SetDirty(orbit);
		}	


	}
}
