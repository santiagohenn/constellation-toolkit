package constellation.tools.math;

import satellite.tools.utils.Utils;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * This class holds mathematical and object transformations
 * **/
public class Transformations {

    private static Map<Double, Double> radii = new HashMap<>();
    private static final boolean USE_CONFORMAL_LATITUDE = false;

    public static List<double[]> toEuclideanPlane(List<double[]> nonEuclideanPolygon, double referenceLat, double referenceLon) {

        radii.clear();

        List<double[]> euclideanPolygon = new ArrayList<>();

        // Transform reference GCS coordinates to radians
        final double referenceLatRads = Math.toRadians(referenceLat);
        final double referenceLonRads = Math.toRadians(referenceLon);

        int segment = 0;

        for (double[] nePair : nonEuclideanPolygon) {

            double lat = Math.toRadians(nePair[0]);
            double lon = Math.toRadians(nePair[1]);
            double localRadius = Utils.EARTH_RADIUS_AVG_KM;
            if (USE_CONFORMAL_LATITUDE) localRadius = computeLocalRadius(lat);

            double k = (2 * localRadius) / (1 + Math.sin(referenceLatRads) * Math.sin(lat) +
                    Math.cos(referenceLatRads) * Math.cos(lat) * Math.cos(lon - referenceLonRads));

            double xStereo = k * Math.cos(lat) * Math.sin(lon - referenceLonRads);
            double yStereo = k * (Math.cos(referenceLatRads) * Math.sin(lat) - Math.sin(referenceLatRads) * Math.cos(lat) * Math.cos(lon - referenceLonRads));

            euclideanPolygon.add(new double[]{xStereo, yStereo});

            radii.put(yStereo, localRadius);

        }

        return euclideanPolygon;

    }

    public static List<double[]> toNonEuclideanPlane(List<double[]> euclideanPolygon, double referenceLat, double referenceLon) {

        List<double[]> GCSPolygon = new ArrayList<>();

        // Transform to radians
        double referenceLatRads = Math.toRadians(referenceLat);
        double referenceLonRads = Math.toRadians(referenceLon);

        for (double[] ePair : euclideanPolygon) {

            double xStereo = ePair[0];
            double yStereo = ePair[1];

            double rho = Math.sqrt(Math.pow(xStereo, 2.000) + Math.pow(yStereo, 2.000));

            double localRadius = Utils.EARTH_RADIUS_AVG_KM;
            if (USE_CONFORMAL_LATITUDE) localRadius = radii.getOrDefault(yStereo, Utils.EARTH_RADIUS_AVG_KM);   // FIXME

            double c = 2 * Math.atan2(rho, 2.0 * localRadius);
            double lat = Math.asin(Math.cos(c) * Math.sin(referenceLatRads) + (yStereo * Math.sin(c) * Math.cos(referenceLatRads)) / rho);
            double lon;

            // For exactly the poles, avoid indeterminate points in the equations
            if (referenceLat == 90.0) {
                lon = referenceLonRads + Math.atan2(xStereo, (-yStereo));
            } else if (referenceLat == -90.0) {
                lon = referenceLonRads + Math.atan2(xStereo, yStereo);
            } else {
                lon = referenceLonRads + Math.atan2((xStereo * Math.sin(c)), (rho * Math.cos(referenceLatRads) * Math.cos(c)
                        - yStereo * Math.sin(referenceLatRads) * Math.sin(c)));
            }

            // Go back to degrees
            lat = Math.toDegrees(lat);
            lon = Math.toDegrees(lon);

            while (lon < -180D) lon += 360;
            while (lon > 180D) lon -= 360;

            GCSPolygon.add(new double[]{lat, lon});

        }

        return GCSPolygon;

    }

    private static double computeLocalRadius(double lat) {

        double latRads = Math.toRadians(lat);
        double confLat = 2 * Math.atan(Math.tan(Math.PI / 4 + latRads / 2) * (Math.pow((1 - Math.E * Math.sin(latRads)) / (1 + Math.E * Math.sin(latRads)), (Math.E / 2)))) - Math.PI / 2;
        return (Utils.EARTH_RADIUS_EQ_M * Math.cos(latRads)) / ((1 - Math.pow(Math.E, 2) * Math.pow(Math.sin(latRads), 2)) * Math.cos(confLat));

    }

    private double conformal2lat(double confLat) {

        return confLat + (Math.pow(Math.E, 2) / 2 + 5 * Math.pow(Math.E, 4) / 24 + 13 * Math.pow(Math.E, 8) / 360) * Math.sin(2 * confLat)
                + (7 * Math.pow(Math.E, 2) / 48 + 29 * Math.pow(Math.E, 4) / 240 + 811 * Math.pow(Math.E, 8) / 11520) * Math.sin(4 * confLat)
                + (7 * Math.pow(Math.E, 6) / 120 + 81 * Math.pow(Math.E, 8) / 1120) * Math.sin(6 * confLat)
                + (4279 * Math.pow(Math.E, 8) / 161280) * Math.sin(8 * confLat);

    }


    public static synchronized List<Pair> doubleList2pairList(List<double[]> polygon) {

        List<Pair> pairList = new ArrayList<>();

        for (double[] coordinate : polygon) {
            pairList.add(new Pair(coordinate[0], coordinate[1]));
        }

        return pairList;

    }

}
