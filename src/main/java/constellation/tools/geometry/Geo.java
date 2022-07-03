package constellation.tools.geometry;

import constellation.tools.math.Pair;
import net.sf.geographiclib.*;
import satellite.tools.utils.Log;
import satellite.tools.utils.Utils;

import java.io.*;
import java.util.ArrayList;
import java.util.List;

/**
 * Contains Geographic operations
 */
public class Geo {

    /**
     * Checks whether the provided FOV contains any of Earth's poles
     *
     * @param FOV       A List of Regions to check
     * @param lambdaMax The Maximum Earth Central Angle of the regions to check
     * @return 0 if the FOV does not contain either the north or South Pole, 1 if it contains the North Pole,
     * -1 if it contains the South Pole
     **/
    public static int checkPoleInclusion(FOV FOV, double lambdaMax) {

        if (computeGeodesic(FOV.getReferenceLat(), FOV.getReferenceLon(), 90, 0) <= lambdaMax) {
            return 1;
        } else if (computeGeodesic(FOV.getReferenceLat(), FOV.getReferenceLon(), -90, 0) <= lambdaMax) {
            return -1;
        }
        return 0;
    }

    /**
     * Computes the angular distance over the euclidean plane given two pair of Region objects. Said objects
     * must have reference coordinates, distance will be computed among those coordinates
     *
     * @param r1 the first FOV
     * @param r2 the second FOV
     * @return Double the computed angular distance in degrees
     **/
    public static double computeGeodesic(FOV r1, FOV r2) {
        return computeGeodesic(r1.getReferenceLat(), r1.getReferenceLon(), r2.getReferenceLat(), r2.getReferenceLon());
    }

    /**
     * Computes the angular distance over the euclidean plane given two pair of Coordinates
     *
     * @param lat1 the first FOV's latitude
     * @param lon1 the first FOV's longitude
     * @param lat2 the second FOV's latitude
     * @param lon2 the second FOV's longitude
     * @return a double value for the computed angular distance, in degrees
     **/
    public static double computeGeodesic(double lat1, double lon1, double lat2, double lon2) {
        GeodesicData g = Geodesic.WGS84.Inverse(lat1, lon1, lat2, lon2,
                GeodesicMask.DISTANCE);
        return g.a12;
    }

    /**
     * This method computes the geodetic area of a list of coordinates, given by Pair objects depicting the polygon's
     * vertices. This method uses the net.sf.geographiclib library.
     *
     * @param pairList A List containing the coordinates of the polygon
     * @return Double the computed area in meters squared
     **/
    public static double computeNonEuclideanSurface(List<Pair> pairList) {

        PolygonArea polygonArea = new PolygonArea(Geodesic.WGS84, false);

        for (Pair pair : pairList) {
            polygonArea.AddPoint(pair.lat, pair.lon);
        }

        PolygonResult result = polygonArea.Compute();
        return Math.abs(result.area);

    }

    /**
     * This method computes the geodetic area of a list of coordinates, given by Pair objects depicting the polygon's
     * vertices. This method uses the net.sf.geographiclib library.
     *
     * @param pairList A List containing the coordinates of the polygon
     * @return Double the computed area in meters squared
     **/
    public static double computeNonEuclideanSurface2(List<double[]> pairList) {

        PolygonArea polygonArea = new PolygonArea(Geodesic.WGS84, false);

        for (double[] pair : pairList) {
            polygonArea.AddPoint(pair[0], pair[1]);
        }

        PolygonResult result = polygonArea.Compute();
        return Math.abs(result.area);

    }

    /**
     * Returns the maximum Lambda for a circular (or otherwise not specified eccentricity) orbit, which is defined as
     * the maximum Earth Central Angle or half of a satellite's "cone FOV" over the surface of the Earth.
     *
     * @param semiMajorAxis       the orbit's semi major axis in meters
     * @param visibilityThreshold the height above horizon visibility threshold in degrees
     * @return Te maximum Earth Central Angle for the access area, in degrees
     **/
    public static double getLambdaMax(double semiMajorAxis, double visibilityThreshold) {
        return getLambdaMax(semiMajorAxis, 0, visibilityThreshold);
    }

    /**
     * This method returns the maximum Lambda, which is defined as the maximum Earth Central Angle or
     * half of a satellite's "cone FOV" over the surface of the Earth.
     *
     * @param semiMajorAxis       the orbit's semi major axis in meters
     * @param eccentricity        the orbit's eccentricity
     * @param visibilityThreshold the height above horizon visibility threshold in degrees
     * @return Te maximum Earth Central Angle for the access area, in degrees
     **/
    public static double getLambdaMax(double semiMajorAxis, double eccentricity, double visibilityThreshold) {

        double hMax = ((1 + eccentricity) * semiMajorAxis) - Utils.EARTH_RADIUS_EQ_M;
        double etaMax = Math.asin((Utils.EARTH_RADIUS_EQ_M * Math.cos(Math.toRadians(visibilityThreshold))) / (Utils.EARTH_RADIUS_EQ_M + hMax));
        return 90 - visibilityThreshold - Math.toDegrees(etaMax);

    }

    /**
     * This method returns a polygon that approximates a circular access area centered at (centerLat, centerLon) with
     * a radius of lambdaMax and the passed number of segments
     *
     * @param lambdaMax the Maximum Earth Central Angle
     * @param centerLat the latitude of the circle's center
     * @param centerLon the longitude of the circle's center
     * @param segments  the amount of segments for the polygon
     * @return a List of double[] containing the polygon (counter clock-wise direction)
     **/
    public static List<double[]> drawCircularAAP(double lambdaMax, double centerLat, double centerLon, double segments) {

        List<double[]> coordinates = new ArrayList<>();

        double lambdaMaxRads = Math.toRadians(lambdaMax);
        double theta;

        // Obtain center latitude' in radians
        double centerLat_ = (Math.PI / 2.0) - (Math.toRadians(centerLat));

        for (int segment = 1; segment <= segments; segment++) {

            // Get the sweep angle
            theta = (segment * 2 * Math.PI) / segments;

            // Obtain current Latitude'
            double pointerLat_ = Math.acos(Math.cos(lambdaMaxRads) * Math.cos(centerLat_)
                    + Math.sin(lambdaMaxRads) * Math.sin(centerLat_) * Math.cos(theta));

            double pointerLat = Math.toDegrees((Math.PI / 2.0) - pointerLat_);

            // Obtain Longitude
            double H = 1;
            double cl_mod_360 = Math.toDegrees(centerLat_) % 360;
            if (180 <= cl_mod_360 && cl_mod_360 < 360) {
                H = -1;
            }

            double deltaLonArgument = ((Math.cos(lambdaMaxRads) - Math.cos(pointerLat_) * Math.cos(centerLat_))
                    / (H * Math.sin(centerLat_) * Math.sin(pointerLat_)));

            // Fix for precision problem near theta ~ PI
            if (deltaLonArgument > 1.0000) {
                deltaLonArgument = 1.0000;
            } else if (deltaLonArgument < -1.0000) {
                deltaLonArgument = -1.00000;
            }

            double deltaLon = Math.acos(deltaLonArgument) - (Math.PI / 2) * (H - 1);

            double pointerLon = centerLon + Math.toDegrees(deltaLon);
            if (theta <= Math.PI) {
                pointerLon = centerLon - Math.toDegrees(deltaLon);
            }

            // Correct special cases
            if (centerLat == -90.0) {
                pointerLon = Math.toDegrees(theta);
                pointerLat = -90.0 + lambdaMax;
            } else if (centerLat == 90.0) {
                pointerLon = Math.toDegrees(theta);
                pointerLat = 90.0 - lambdaMax;
            }

            if (theta == 2 * Math.PI && (centerLat + lambdaMax) == 90.0) {
                pointerLat = 90.0;
                pointerLon = 0;
            } else if (theta == 2 * Math.PI && (centerLat + lambdaMax) == -90.0) {
                pointerLat = -90.0;
                pointerLon = 0;
            }

            while (pointerLon < -180D) {
                pointerLon = pointerLon + 360;
            }

            while (pointerLon > 180D) {
                pointerLon = pointerLon - 360;
            }

            if (Double.isNaN(pointerLat) || Double.isNaN(pointerLon)) {

                Log.debug("Segment: " + segment + " Center coordinates: (" + centerLat + "," + centerLon + ")");
                Log.debug("centerLat_ = " + centerLat_ + " - centerLatDeg_ = " + Math.toDegrees(centerLat_));
                Log.debug("Math.toDegrees(centerLat_) % 360 = " + (Math.toDegrees(centerLat_) % 360));
                Log.debug("acos argument: " + deltaLonArgument);
                Log.debug("acos: " + Math.acos(deltaLonArgument));
                Log.debug("deltaLon = " + deltaLon);
                Log.debug("pointerLon = " + pointerLon);
                Log.debug("theta: " + theta);
                Log.debug("Math.cos(lambdaMaxRads) = " + Math.cos(lambdaMaxRads)
                        + " Math.cos(pointerLat_) = " + Math.cos(pointerLat_) + " Math.cos(centerLat_) = " + Math.cos(centerLat_)
                        + " Math.sin(pointerLat_) = " + Math.sin(pointerLat_));

            }

            coordinates.add(new double[]{pointerLat, pointerLon});

        }

        return coordinates;

    }

    /**
     * Reads a file containing asset(s) parameter(s) and returns a list of objects accordingly
     *
     * @return List<Pair>
     */
    public static List<double[]> file2DoubleList(String fileName) {

        List<double[]> pairList = new ArrayList<>();
        var file = new File(fileName);
        try (var fr = new FileReader(file); var br = new BufferedReader(fr)) {
            String line;
            while ((line = br.readLine()) != null) {
                if (!line.startsWith("//") && line.length() > 0) {
                    var data = line.split(",");
                    pairList.add(new double[]{Double.parseDouble(data[0]), Double.parseDouble(data[1])});
                }
            }
        } catch (FileNotFoundException e) {
            Log.error("Unable to find file: " + fileName);
            e.printStackTrace();
        } catch (IOException e) {
            Log.error("IOException: " + fileName);
            e.printStackTrace();
        }

        return pairList;
    }


}
