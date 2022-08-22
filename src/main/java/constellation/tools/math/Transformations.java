package constellation.tools.math;

import java.util.ArrayList;
import java.util.List;

import static java.lang.Math.PI;

/**
 * This class holds mathematical and object transformations
 **/
public class Transformations {

    public static final double EARTH_RADIUS_EQ_M = 6378135.0D;
    private static final double RADIUS = 6371.009D;
    public static final double WGS84_F = 1 / 298.257223563;
    public static final double WGS84_E2 = 0.00669437999014;
    public static final double WGS84_E = 0.081819190842613;

    public static List<double[]> toEuclideanPlane(List<double[]> nonEuclideanPolygon, double referenceLat, double referenceLon) {

        List<double[]> euclideanPolygon = new ArrayList<>();

        // Transform reference GCS coordinates to radians
        final double referenceLatRads = Math.toRadians(referenceLat);
        final double referenceLonRads = Math.toRadians(referenceLon);

        for (double[] nePair : nonEuclideanPolygon) {

            double lat = Math.toRadians(nePair[0]);
            double lon = Math.toRadians(nePair[1]);

            euclideanPolygon.add(toStereo(lat, lon, referenceLatRads, referenceLonRads));

        }

        return euclideanPolygon;

    }

    /**
     * Takes the a pair of geographic coordinates and transforms them into the euclidean plane
     **/
    public static double[] toStereo(double lat, double lon, double referenceLatRads, double referenceLonRads) {

        double localR = RADIUS; // getLocalRadiusKm(Math.toRadians(lat));
        double k = (2 * localR) / (1 + Math.sin(referenceLatRads) * Math.sin(lat) +
                Math.cos(referenceLatRads) * Math.cos(lat) * Math.cos(lon - referenceLonRads));

        double xStereo = k * Math.cos(lat) * Math.sin(lon - referenceLonRads);
        double yStereo = k * (Math.cos(referenceLatRads) * Math.sin(lat) - Math.sin(referenceLatRads) * Math.cos(lat) * Math.cos(lon - referenceLonRads));

        return new double[]{xStereo, yStereo};

    }

    public static List<double[]> toNonEuclideanPlane(List<double[]> euclideanPolygon, double referenceLat, double referenceLon) {

        List<double[]> GCSPolygon = new ArrayList<>();

        // Transform to radians
        double referenceLatRads = Math.toRadians(referenceLat);
        double referenceLonRads = Math.toRadians(referenceLon);

        for (double[] ePair : euclideanPolygon) {

            double[] nePair = toGeodetic(ePair[0], ePair[1], referenceLatRads, referenceLonRads);
            GCSPolygon.add(new double[]{nePair[0], nePair[1]});

        }

        return GCSPolygon;

    }

    public static double[] toGeodetic(double xStereo, double yStereo, double referenceLatRads, double referenceLonRads) {

        double rho = Math.sqrt(Math.pow(xStereo, 2.000) + Math.pow(yStereo, 2.000));

        double c = 2 * Math.atan2(rho, 2.0 * RADIUS);
        double lat = Math.asin(Math.cos(c) * Math.sin(referenceLatRads) + (yStereo * Math.sin(c) * Math.cos(referenceLatRads)) / rho);
        double lon;

        // For exactly the poles, avoid indeterminate points in the equations
        if (referenceLatRads == PI / 2) {
            lon = referenceLonRads + Math.atan2(xStereo, (-yStereo));
        } else if (referenceLatRads == -PI / 2) {
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

        return new double[]{lat, lon};

    }

    public static double getLocalRadiusKm(double sphericalLat) {
        double cLat = lat2conformal(sphericalLat);
        return ((EARTH_RADIUS_EQ_M * Math.cos(Math.toRadians(sphericalLat)))
                / ((1 - WGS84_E2 * Math.pow(sphericalLat, 2)) * Math.cos(cLat))) / 1000.0;
    }

    public static double lat2ConformalD(double confLat) {
        return Math.toDegrees(lat2conformal(Math.toRadians(confLat)));
    }

    public static double lat2conformal(double sphericalLat) {
        return sphericalLat + (Math.pow(WGS84_E, 2) / 2 + 5 * Math.pow(WGS84_E, 4) / 24 + 3 * Math.pow(WGS84_E, 6) / 32 + 281 * Math.pow(WGS84_E, 8) / 5760) * Math.sin(2 * sphericalLat)
                + (5 * Math.pow(WGS84_E, 4) / 48 + 7 * Math.pow(WGS84_E, 6) / 80 + 697 * Math.pow(WGS84_E, 8) / 11520) * Math.sin(4 * sphericalLat)
                + (13 * Math.pow(WGS84_E, 6) / 480 + 461 * Math.pow(WGS84_E, 8) / 13440) * Math.sin(6 * sphericalLat)
                + (1237 * Math.pow(WGS84_E, 8) / 161280) * Math.sin(8 * sphericalLat);
    }

    public double conformal2lat(double confLat) {

        return confLat + (Math.pow(WGS84_E, 2) / 2 + 5 * Math.pow(WGS84_E, 4) / 24 + 1 * Math.pow(WGS84_E, 6) / 12 + 13 * Math.pow(WGS84_E, 8) / 360) * Math.sin(2 * confLat)
                + (7 * Math.pow(WGS84_E, 4) / 48 + 29 * Math.pow(WGS84_E, 6) / 240 + 811 * Math.pow(WGS84_E, 8) / 11520) * Math.sin(4 * confLat)
                + (7 * Math.pow(WGS84_E, 6) / 120 + 81 * Math.pow(WGS84_E, 8) / 1120) * Math.sin(6 * confLat)
                + (4279 * Math.pow(WGS84_E, 8) / 161280) * Math.sin(8 * confLat);

    }

    public static double unix2julian(long unix) {
        return ( unix / 86400.0 ) + 2440587.5;
    }

}
