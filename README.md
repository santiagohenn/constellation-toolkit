# constellation-toolkit
A set of tools to analyze coverage areas served by multiple Low-Earth Orbit satellites.

## Dependencies

You need to have the Orekit dataset installed on your computer. You can clone it from [HERE](https://github.com/santiagohenn/orekit-data.git) or find the latest version on the [official site](https://orekit.org).

## Basic usage

The easiest approach is to modify the configuration file with the scenario settings and run the app:

```java
String configurationsPath = "C:/experiments/config.properties";
ConstellationCoverageComputer app = new ConstellationCoverageComputer(configurationsPath);
app.run();
```

An example properties file with explanations is provided in this repo. Here are some formats for the required inputs.

#### Satellite constellation

A CSV file containing a satellite on each line, with no headers. Each line contains the following information, in this order:

```tex
Elements' date (ISO 8601), Semi-major axis (m), Eccentricity, Inclination (°), RAAN (°), Arg. of Periapsis (°), True Anomaly (°)
```
_Example:_

```csv
2025-01-01T16:00:00.000,6978135,0,45,250,0,220
2025-01-01T16:00:00.000,6978135,0,45,250,0,240
2025-01-01T16:00:00.000,6978135,0,45,250,0,260
2025-01-01T16:00:00.000,6978135,0,45,250,0,280
```

#### Polygon (ROI) File

A CSV file containing coordinates that define the Region of Interest (ROI). Each line contains coordinates in the format `latitude,longitude`, in degrees.

_Example:_

```csv
-24.492,-67.165
-24.492,-66.743
-24.492,-66.322
...
```

Important: Do not use self-intersecting polygons. Use clockwise order for the coordinates according to the surface you want to enclose. Points will be placed at surface height according to the WGS84 Earth model. Huge polygons (e.g., a complete Earth hemisphere) might yield strange results.

## Output format

The tool outputs several files. The coverage output file is a time series where the first column is the epoch time in milliseconds and the next columns are the percentage of ROI accessed by at least `1,2,...,N` satellites.

_example:_

```csv
...
900000,0.0,0.0,0.0
930000,0.0,0.0,0.0
960000,2.268035445461762,0.0,0.0
990000,8.93095638620536,0.0,0.0
1020000,17.758691643704537,0.0,0.0
1050000,27.99764177061209,0.0,0.0
1080000,39.17666740891235,3.9883331407036344,0.0
1110000,50.902285625179864,11.417244533935357,0.0
1140000,62.79103986849493,20.718427471311266,0.0
1170000,74.41302373179948,31.275092937720313,0.575558227330941
1200000,85.22351867972844,42.64554395335318,5.979435751486388
1230000,94.38078799473708,54.4444032137197,14.04790437901466
1260000,99.96486179896316,66.2822145849954,23.771814105655185
1290000,100.0,77.70474322366347,34.601016370118366
1320000,100.0,88.120667225942,46.11918388728884
1350000,100.0,96.49243601128394,57.94926750115989
1380000,100.0,100.0,69.6926009314914
1410000,100.0,100.0,80.86883278458896
1440000,100.0,100.0,90.81400701447939
1470000,100.0,100.0,98.26693778785267
1500000,100.0,100.0,100.0
1530000,100.0,100.0,100.0
1560000,100.0,100.0,100.0
1590000,100.0,100.0,100.0
...
```

Coverage results are also given in JSON format, both for the entire constellation and for each satellite. If geographic outputs are enabled, more JSON files are generated with coordinates over time and for a particular snapshot time, provided this is enabled in the configuration file.

### Constellation hash

Each constellation is assigned a unique hash based on the orbital elements of its satellites. This hash is included in the output files for identification purposes. The hashing process works as follows:

- For each satellite, all orbital elements except the semi-major axis are multiplied by 1000 (the semi-major axis remains unchanged).
- The integer part of each element is taken.
- The elements are ordered as: semi-major axis, inclination, eccentricity, RAAN, argument of periapsis, and true anomaly.
- The integer values are concatenated into a single string for each satellite, with no separators.
- All satellite strings are sorted in natural (lexicographical) order and concatenated into one long string.
- This final string is hashed using the SHA-256 algorithm to produce the constellation hash.

An implementation for hashing constellation parameters in Python is included in the helper_scripts folder.

## Advanced usage

Most of the configurations can be overridden at runtime. This is useful for scenario optimization, parameter exploration, etc.

Change the output path:

```java
app.setOutputPath(outputPath);
```

Change the constellation from a file:

```java
List<Satellite> satelliteList = Utils.satellitesFromFile(ConstellationFilePath);
app.setSatelliteList(satelliteList);
```
Change a single satellite element(s):

```java
List<Satellite> satelliteList = app.getSatelliteList();
satelliteList.get(0).getElements().setInclination(98.888);
satelliteList.get(0).getElements().setSemiMajorAxis(6771000);
```
_The satellite index corresponds to its zero-based position in the constellation file._

Change the ROI:

```java
List<double[]> ROI = FileUtils.file2DoubleList("polygon_path");
app.setROI(ROI);
```

You can generate a polygon that encloses a spherical cap on the Earth and set it as the scenario ROI:

```java
List<double[]> ROI = GeographicTools.computeSphericalCap(centerLatDegrees, centerLonDegrees, roiSphericalRadius, roiSegments);
app.setROI(ROI);
```

If you need to compute some geographic parameters, `GeographicTools` has some useful methods. For example, to compute the surface area (in square meters):

```java
GeographicTools geographicTools = new GeographicTools();
double roiSurface = geographicTools.computeNonEuclideanSurface(ROI);
```

### Properties detailed explanation

```properties
# Path to the output folder (e.g., C:/Projects/coverage-analysis/output)
output_path=C:/Projects/coverage-analysis/output

# Path to Orekit data (e.g., C:/Orekit/data)
orekit_data_path=C:/orekit-data

# Path to satellite orbital elements file (CSV format)
satellites_file=C:/Projects/constellations/satellites.csv

# Whether to use a ROI (Region of Interest) polygon file (true or false)
use_roi_file=true

# Path to the ROI polygon file
roi_path=C:/Projects/rois/argentina.csv

# Scenario time parameters

# Start date and time in ISO 8601 format (e.g., 2025-01-01T16:00:00.000)
start_date=2025-01-01T16:00:00.000

# End date and time in ISO 8601 format (e.g., 2025-01-01T17:30:00.000)
end_date=2025-01-01T17:30:00.000

# Time step in seconds (e.g., 30)
time_step_seconds=30

# Whether to propagate internally (true) or use external positions (false)
propagate_internally=true

# Path to the directory with the satellite ephemeris (only if propagate_internally=false)
positions_path=C:/Projects/constellation-toolkit/positions

# Visibility threshold in degrees (e.g., 20)
visibility_threshold=20

# Maximum number of interacting satellites at any time (e.g., 3)
max_subset_size=3

# Number of segments to approximate the polygon for the ROI (e.g., 50 for a good approximation)
polygon_segments=50

# Whether to save geographic data in the output (true or false)
# This is useful for debugging and visualization, but can increase output size significantly
save_geographic=false

# Takes a "picture" at milliseconds from epoch of all the satellite's FOVs and intersections (e.g., 86400) leave blank for no snapshot
snapshot_time_milliseconds=

#####################################################################################
# Advanced configs (try to keep these as they are unless you know what you are doing)

# Enable debug mode (true or false). This will log A LOT of data.
debug_mode=false
# Threshold detection value (e.g., 0.001). Used for the polygon intersector.
th_detection=0.001
# Polygon epsilon value (e.g., 0.0001). Used for the polygon intersector.
polygon_epsilon=0.0001
# Lambda exclusion value (e.g., 3). Used for performance improvement.
lambda_exclusion=3
# Whether to save Euclidean data (true or false), like geo data but in the euclidean plane.
save_euclidean=true
#####################################################################################
```