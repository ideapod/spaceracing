On-line documentation/tutorial videos:
http://nbodyphysics.com/blog/gravity-engine-doc-1-3-2/

Docs for script elements:
http://nbodyphysics.com/gravityengine/html/

Support: nbodyphysics@gmail.com

Gravity Engine 2.2 (Nov 6, 2018)
================================

Features:
1) Add code for Lambert transfer to a point. 

2) Add a ManeuverRenderer to see velocity changes for Lambert maneuvers. 

3) Add OrbitSegment renderer to show which way a Lambert transfer will go (short/long)

4) Add GetMass() function to retrieve internal physics engine mass. 

5) Convert PatchedConicXfer internals to double

Bugs:
1) Fix FixedObject to take position from transform when units = DIMENSIONLESS

2) Position/velocity precalc when parent object has a Kepler OrbitEllipse when application is playing. 

3) Implement deltaV for LambertMinEnergy and LambertUniversal.

4) Remove LambertMinEnergy, since LambertUniversal does this as first pass and this removes duplicate code. 
(LambertMinEnergy had a bug and it did not make sense to fix it.)

5) Issue #23 (Exception on Kepler object on start). Breakage due to (2) but also likely existed in some
cases before this. Depended on order of object insertion. 

6) Correct gizmo display of hill sphere in editor view when in scaled units. 
[Issue #24]

7) Fix loop convergence check in LambertUniversal code.


Gravity Engine 2.1 Sept. 11,2018
================================

New:
1) Add support for a selective (per body) custom force.

2) Add a full solar system scene. (My solar system includes Pluto!)


Bug Fixes:
1) Allow editor float field drag by removing DelayedFloatField from: 
	EllipseBaseEditor
	NBodyEditor
	RandomPlanetEditor
	SolarBodyEditor
	SolarSystemEditor

2) Fix delete of binary planets and moons in AddDeleteTester. Add a Chaos Monkey checkbox for random testing. 

3) Fix adjustment to eccentricity in SolarSystemBody (could not be changed in the editor)

4) Modify LambertDemoController to ensure key presses are limited to scene states that are relevent. 

5) Fix outer to inner Lambert transfer. Sone issues remain

6) Adjust parameters in the Tutorial 6. 
- reduce thrust
- uncheck align motion with trajectory

7) Fix collision detection bug. 

Known Issues:
=============
1) Lambert Xfers
Outer to inner transfers sometimes fail to match the inner orbit
In some cases the selected direction on the transfer ellipse is not optimal (wrong direction)

2) AddDeleteTester
When the "Run Chaos Monkey" is checked get sporadic AABB and isFinite errors in short bursts. 
This is related to the OrbitPredictors and is under investigation.  


Gravity Engine 2.0
August 12, 2018
================

Lambert Orbital Transfers 
Lunar Course Corrections
GEConsole command additions
Add massless option to FrameRate tester
Support Gravity State evolution on an independent thread

Fixes:
- hyperbola course predition fixes
- ordering of updates for nested Kepler objects



On-line documentation/tutorial videos:
http://nbodyphysics.com/blog/gravity-engine-doc-1-3-2/

Docs for script elements:
http://nbodyphysics.com/gravityengine/html/

Support: nbodyphysics@gmail.com

