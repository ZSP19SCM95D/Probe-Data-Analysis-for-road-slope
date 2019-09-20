1. Partition6467ProbePoints.csv and Partition6467LinkData.csv are two Input files.  Partition6467MatchedPoints.csv is the mapmatching files.
Prtition6467Slopes.csv is the Slope data we get


2. ProbePoints Record Format:

	sampleID, dateTime, sourceCode, latitude, longitude, altitude, speed, heading


   LinkData Record Format:

	linkPVID, refNodeID, nrefNodeID, length, functionalClass, directionOfTravel, speedCategory, fromRefSpeedLimit, toRefSpeedLimit, fromRefNumLanes, toRefNumLanes, multiDigitized, urban, timeZone, shapeInfo, curvatureInfo, slopeInfo(in decimal degrees) and elevation (in decimal meters)


MatchedPoints Record Format:

	sampleID, dateTime, sourceCode, latitude, longitude, altitude, speed, heading, linkPVID, direction, distFromRef, distFromLink

3. The way to run the program:

Python Probe.py 

Make sure that the two input files are in the same directory of Probe.py


 