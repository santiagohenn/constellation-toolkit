package constellation.tools.geometry;

import net.sf.geographiclib.*;
import org.apache.commons.math3.geometry.euclidean.threed.Vector3D;
import satellite.tools.utils.Utils;

import java.util.ArrayList;
import java.util.List;

/**
 * Contains Geographic operations
 */
public class GeographicTools {

    public static final double EARTH_RADIUS_EQ_M = 6378135.0D;
    public static final double WGS84_EQ_RADIUS_M = 6378137.0D;
    public static final double WGS84_EQ_RADIUS_KM = 6378.137D;
    public static final double WGS84_F = 1 / 298.257223563;
    public static final double WGS84_E2 = 0.00669437999014;
    public static final double WGS84_E = 0.081819190842613;

    private static final long POLAR_NO_ERROR = 0x0000;
    private static final long POLAR_LAT_ERROR = 0x0001;
    private static final long POLAR_LON_ERROR = 0x0002;
    private static final long POLAR_ORIGIN_LAT_ERROR = 0x0004;
    private static final long POLAR_ORIGIN_LON_ERROR = 0x0008;
    public static final long POLAR_EASTING_ERROR = 0x0010;
    public static final long POLAR_NORTHING_ERROR = 0x0020;
    private static final long POLAR_A_ERROR = 0x0040;
    private static final long POLAR_INV_F_ERROR = 0x0080;
    public static final long POLAR_RADIUS_ERROR = 0x0100;

    private static final double PI = 3.14159265358979323;
    private static final double PI_OVER_2 = PI / 2.0;
    private static final double PI_OVER_4 = PI / 4.0;
    private static final double TWO_PI = 2.0 * PI;

    /* Ellipsoid Parameters, default to WGS 84  */
    private double polarA = 6378137.0;                    /* Semi-major axis of ellipsoid in meters  */
    private double polarF = 1 / 298.257223563;            /* Flattening of ellipsoid  */
    private double es = 0.08181919084262188000;           /* Eccentricity of ellipsoid    */
    private double esOver2 = .040909595421311;            /* es / 2.0 */
    private double southernHemisphere = 0;                /* Flag variable */
    private double mc = 1.0;
    private double tc = 1.0;
    private double e4 = 1.0033565552493;
    private double polarAMc = 6378137.0;                 /* Polar_a * mc */
    private double twoPolarA = 12756274.0;               /* 2.0 * Polar_a */

    /* Polar Stereographic projection Parameters */
    private double polarOriginLat = ((PI * 90) / 180);   /* Latitude of origin in radians */
    private double polarOriginLong = 0.0;                /* Longitude of origin in radians */
    private double polarFalseEasting = 0.0;              /* False easting in meters */
    private double polarFalseNorthing = 0.0;             /* False northing in meters */

    /* Maximum variance for easting and northing values for WGS 84. */
    private double polarDeltaEasting = 12713601.0;
    private double polarDeltaNorthing = 12713601.0;

    private double easting;
    private double northing;
    private double latitude;
    private double longitude;

    /**
     * Checks whether the provided FOV contains any of Earth's poles
     *
     * @param accessRegion A List of Regions to check
     * @param lambdaMax    The Maximum Earth Central Angle of the regions to check
     * @return 0 if the FOV does not contain either the north or South Pole, 1 if it contains the North Pole,
     * -1 if it contains the South Pole
     **/
    public int checkPoleInclusion(AccessRegion accessRegion, double lambdaMax) {
        if (computeGeodesic(accessRegion.getSspLat(), accessRegion.getSspLon(), 90, 0) <= lambdaMax) {
            return 1;
        } else if (computeGeodesic(accessRegion.getSspLat(), accessRegion.getSspLon(), -90, 0) <= lambdaMax) {
            return -1;
        }
        return 0;
    }

    /**
     * Checks whether the provided refLat and refLon could contain any of Earth's poles
     *
     * @param refLat    The reference
     * @param refLon    A List of Regions to check
     * @param lambdaMax The Maximum Earth Central Angle of the regions to check
     * @return 0 if the FOV does not contain either the north or South Pole, 1 if it contains the North Pole,
     * -1 if it contains the South Pole
     **/
    public int checkPoleInclusion(double refLat, double refLon, double lambdaMax) {

        if (computeGeodesic(refLat, refLon, 90, 0) <= lambdaMax) {
            return 1;
        } else if (computeGeodesic(refLat, refLon, -90, 0) <= lambdaMax) {
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
    public double computeGeodesic(AccessRegion r1, AccessRegion r2) {
        return computeGeodesic(r1.getSspLat(), r1.getSspLon(), r2.getSspLat(), r2.getSspLon());
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
    public double computeGeodesic(double lat1, double lon1, double lat2, double lon2) {
        GeodesicData g = Geodesic.WGS84.Inverse(lat1, lon1, lat2, lon2,
                GeodesicMask.DISTANCE);
        return g.a12;
    }

    /**
     * This method computes the geodetic area of a list of coordinates, given by Pair objects depicting the polygon's
     * vertices. This method uses the net.sf.geographiclib library. Paper: https://doi.org/10.1007/s00190-012-0578-z
     * Section 6 , C. F. F. Karney, Algorithms for geodesics, J. Geodesy 87, 43–55 (2013).
     *
     * @param pairList A List containing the coordinates of the polygon
     * @return Double the computed area in meters squared
     **/
    public double computeNonEuclideanSurface(List<double[]> pairList) {

        if (pairList.isEmpty()) {
            return 0;
        }

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
    public double getLambdaMax(double semiMajorAxis, double visibilityThreshold) {
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
    public double getLambdaMax(double semiMajorAxis, double eccentricity, double visibilityThreshold) {

        double hMax = ((1 + eccentricity) * semiMajorAxis) - Utils.EARTH_RADIUS_EQ_M;
        double etaMax = Math.asin((Utils.EARTH_RADIUS_EQ_M * Math.cos(Math.toRadians(visibilityThreshold))) / (Utils.EARTH_RADIUS_EQ_M + hMax));
        return 90 - visibilityThreshold - Math.toDegrees(etaMax);

    }

    /**
     * This method returns the maximum Lambda, which is defined as the maximum Earth Central Angle or
     * half of a satellite's "cone FOV" over the surface of the Earth.
     *
     * @param x                   x component of the position vector in ECEF coordinates
     * @param y                   y component of the position vector in ECEF coordinates
     * @param z                   z component of the position vector in ECEF coordinates
     * @param visibilityThreshold the height above horizon visibility threshold in degrees
     * @return Te maximum Earth Central Angle for the access area, in degrees
     **/
    public double getLambdaMax(double x, double y, double z, double visibilityThreshold) {

        Vector3D pos = new Vector3D(x, y, z);
        double hMax = pos.getNorm() - WGS84_EQ_RADIUS_KM;
        double etaMax = Math.asin((WGS84_EQ_RADIUS_KM * Math.cos(Math.toRadians(visibilityThreshold))) / (WGS84_EQ_RADIUS_KM + hMax));
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
    public static List<double[]> computeSphericalCap(double lambdaMax, double centerLat, double centerLon, double segments) {

        List<double[]> coordinates = new ArrayList<>();

        double lambdaMaxRads = Math.toRadians(lambdaMax);
        double theta = 0;
        double centerLatPrime = (Math.PI / 2.0) - (Math.toRadians(centerLat));
        double angleStep = ( 2.0 * Math.PI) / (segments - 1);

        for (int segment = 1; segment <= segments; segment++) {

            theta += angleStep;

            // Obtain current Latitude'
            double pointerLatPrime = Math.acos(Math.cos(lambdaMaxRads) * Math.cos(centerLatPrime)
                    + Math.sin(lambdaMaxRads) * Math.sin(centerLatPrime) * Math.cos(theta));

            double pointerLat = Math.toDegrees((Math.PI / 2.0) - pointerLatPrime);

            // Obtain Longitude
            double clMod360 = Math.toDegrees(centerLatPrime) % 360;
            int hParameter = (180 <= clMod360 && clMod360 < 360) ? -1 : 1;

            double deltaLonArgument = ((Math.cos(lambdaMaxRads) - Math.cos(pointerLatPrime) * Math.cos(centerLatPrime))
                    / (hParameter * Math.sin(centerLatPrime) * Math.sin(pointerLatPrime)));

            // Fix for precision problem near theta ~ PI
            deltaLonArgument = Math.min(deltaLonArgument, 1.00000);
            deltaLonArgument = Math.max(deltaLonArgument, -1.00000);

            double deltaLon = Math.acos(deltaLonArgument) - (Math.PI / 2.0) * (hParameter - 1);
            double pointerLon = (theta <= Math.PI) ? centerLon - Math.toDegrees(deltaLon) : centerLon + Math.toDegrees(deltaLon);

            // Correct special cases
            if (centerLat == -90.0) {
                pointerLon = Math.toDegrees(theta);
                pointerLat = -90.0 + lambdaMax;
            } else if (centerLat == 90.0) {
                pointerLon = Math.toDegrees(theta);
                pointerLat = 90.0 - lambdaMax;
            }

            if (theta == 2.0 * Math.PI && (centerLat + lambdaMax) == 90.0) {
                pointerLat = 90.0;
                pointerLon = 0.0;
            } else if (theta == 2.0 * Math.PI && (centerLat + lambdaMax) == -90.0) {
                pointerLat = -90.0;
                pointerLon = 0.0;
            }

            pointerLon = net.sf.geographiclib.GeoMath.AngNormalize(pointerLon);

            coordinates.add(new double[]{pointerLat, pointerLon});

        }
        return coordinates;
    }

    public double[] computeAntipode(double lat, double lon) {

        double[] antipode = new double[]{lat, lon};
        antipode[0] = -1 * antipode[0];
        antipode[1] = 180 - antipode[1];

        while (antipode[1] < -360) {
            antipode[1] += 360;
        }

        while (antipode[1] > 360) {
            antipode[1] -= 360;
        }

        return antipode;

    }

    public List<double[]> toEuclideanPlane(List<double[]> nonEuclideanPolygon, double referenceLat, double referenceLon) {

        List<double[]> euclideanPolygon = new ArrayList<>();

        // Transform reference GCS coordinates to radians
        final double referenceLatRads = Math.toRadians(referenceLat);
        final double referenceLonRads = Math.toRadians(referenceLon);

        setPolarStereographicParameters(WGS84_EQ_RADIUS_KM, WGS84_F, referenceLatRads,
                referenceLonRads, 0, 0);

        for (double[] nePair : nonEuclideanPolygon) {

            double lat = Math.toRadians(nePair[0]);
            double lon = Math.toRadians(nePair[1]);
            convertGeodeticToPolarStereographic(lat, lon);
            euclideanPolygon.add(new double[]{this.easting, this.northing});

        }

        return euclideanPolygon;

    }

    public List<double[]> toNonEuclideanPlane(List<double[]> euclideanPolygon, double referenceLat, double referenceLon) {

        List<double[]> GCSPolygon = new ArrayList<>();

        // Transform to radians
        double referenceLatRads = Math.toRadians(referenceLat);
        double referenceLonRads = Math.toRadians(referenceLon);

        setPolarStereographicParameters(WGS84_EQ_RADIUS_KM, WGS84_F, referenceLatRads,
                referenceLonRads, 0, 0);

        for (double[] ePair : euclideanPolygon) {

            double xStereo = ePair[0];
            double yStereo = ePair[1];

            convertPolarStereographicToGeodetic(xStereo, yStereo);

            GCSPolygon.add(new double[]{Math.toDegrees(this.latitude), Math.toDegrees(this.longitude)});

        }

        return GCSPolygon;

    }

    /**
     * Takes the a pair of geographic coordinates and transforms them into the euclidean plane
     **/
    @Deprecated
    public static double[] toStereo(double lat, double lon, double referenceLatRads, double referenceLonRads) {

        double localRadius = Utils.EARTH_RADIUS_AVG_KM;

        double k = (2 * localRadius) / (1 + Math.sin(referenceLatRads) * Math.sin(lat) +
                Math.cos(referenceLatRads) * Math.cos(lat) * Math.cos(lon - referenceLonRads));

        double xStereo = k * Math.cos(lat) * Math.sin(lon - referenceLonRads);
        double yStereo = k * (Math.cos(referenceLatRads) * Math.sin(lat) - Math.sin(referenceLatRads) * Math.cos(lat) * Math.cos(lon - referenceLonRads));


        return new double[]{xStereo, yStereo};

    }

    /**
     * Takes a pair of geographic coordinates and transforms them into the euclidean plane
     **/
    @Deprecated
    public static double[] toGeodetic(double xStereo, double yStereo, double referenceLatRads, double referenceLonRads) {

        double rho = Math.sqrt(Math.pow(xStereo, 2.000) + Math.pow(yStereo, 2.000));

        double localRadius = Utils.EARTH_RADIUS_AVG_KM;

        double c = 2 * Math.atan2(rho, 2.0 * localRadius);
        double lat = Math.asin(Math.cos(c) * Math.sin(referenceLatRads) + (yStereo * Math.sin(c) * Math.cos(referenceLatRads)) / rho);
        double lon;

        // For exactly the poles, avoid indeterminate points in the equations
        if (referenceLatRads == Math.PI) {
            lon = referenceLonRads + Math.atan2(xStereo, (-yStereo));
        } else if (referenceLonRads == -Math.PI) {
            lon = referenceLonRads + Math.atan2(xStereo, yStereo);
        } else {
            lon = referenceLonRads + Math.atan2((xStereo * Math.sin(c)), (rho * Math.cos(referenceLatRads) * Math.cos(c)
                    - yStereo * Math.sin(referenceLatRads) * Math.sin(c)));
        }

        // Go back to degrees
        lat = Math.toDegrees(lat);
        lon = net.sf.geographiclib.GeoMath.AngNormalize(lon);

        return new double[]{lat, lon};

    }

    /**
     * The function setPolarStereographicParameters receives the ellipsoid parameters and Polar Stereograpic projection
     * parameters as inputs, and sets the corresponding state variables.  If any errors occur, error code(s) are
     * returned by the function, otherwise POLAR_NO_ERROR is returned.
     *
     * @param a                        Semi-major axis of ellipsoid, in meters
     * @param f                        Flattening of ellipsoid
     * @param latitudeOfTrueScale   Latitude of true scale, in radians
     * @param longitudeDownFromPole Longitude down from pole, in radians
     * @param falseEasting            Easting (X) at center of projection, in meters
     * @param falseNorthing           Northing (Y) at center of projection, in meters
     * @return error code
     */
    @SuppressWarnings("squid:S3776")
    public long setPolarStereographicParameters(double a, double f, double latitudeOfTrueScale,
                                                double longitudeDownFromPole, double falseEasting, double falseNorthing) {
        double es2;
        double sLat;
        double cLat;
        double esSin;
        double onePlusEs;
        double oneMinusEs;
        double powEs;
        double invF = 1 / f;
        final double epsilon = 1.0e-2;
        long errorCode = POLAR_NO_ERROR;

        if (a <= 0.0) { /* Semi-major axis must be greater than zero */
            errorCode |= POLAR_A_ERROR;
        }
        if ((invF < 250) || (invF > 350)) { /* Inverse flattening must be between 250 and 350 */
            errorCode |= POLAR_INV_F_ERROR;
        }
        if ((latitudeOfTrueScale < -PI_OVER_2) || (latitudeOfTrueScale > PI_OVER_2)) { /* Origin Latitude out of range */
            errorCode |= POLAR_ORIGIN_LAT_ERROR;
        }
        if ((longitudeDownFromPole < -PI) || (longitudeDownFromPole > TWO_PI)) { /* Origin Longitude out of range */
            errorCode |= POLAR_ORIGIN_LON_ERROR;
        }

        if (errorCode == POLAR_NO_ERROR) { /* no errors */

            polarA = a;
            twoPolarA = 2.0 * polarA;
            polarF = f;

            if (longitudeDownFromPole > PI)
                longitudeDownFromPole -= TWO_PI;
            if (latitudeOfTrueScale < 0) {
                southernHemisphere = 1;
                polarOriginLat = -latitudeOfTrueScale;
                polarOriginLong = -longitudeDownFromPole;
            } else {
                southernHemisphere = 0;
                polarOriginLat = latitudeOfTrueScale;
                polarOriginLong = longitudeDownFromPole;
            }
            polarFalseEasting = falseEasting;
            polarFalseNorthing = falseNorthing;

            es2 = 2 * polarF - polarF * polarF;
            es = Math.sqrt(es2);
            esOver2 = es / 2.0;

            if (Math.abs(Math.abs(polarOriginLat) - PI_OVER_2) > 1.0e-10) {
                sLat = Math.sin(polarOriginLat);
                esSin = es * sLat;
                powEs = Math.pow((1.0 - esSin) / (1.0 + esSin), esOver2);
                cLat = Math.cos(polarOriginLat);
                mc = cLat / Math.sqrt(1.0 - esSin * esSin);
                polarAMc = polarA * mc;
                tc = Math.tan(PI_OVER_4 - polarOriginLat / 2.0) / powEs;
            } else {
                onePlusEs = 1.0 + es;
                oneMinusEs = 1.0 - es;
                e4 = Math.sqrt(Math.pow(onePlusEs, onePlusEs) * Math.pow(oneMinusEs, oneMinusEs));
            }
        }

        /* Calculate Radius */
        convertGeodeticToPolarStereographic(0, polarOriginLong);

        polarDeltaNorthing = northing * 2; // Increased range for accepted easting and northing values
        polarDeltaNorthing = Math.abs(polarDeltaNorthing) + epsilon;
        polarDeltaEasting = polarDeltaNorthing;

        return (errorCode);
    }

    /**
     * The function Convert_Geodetic_To_Polar_Stereographic converts geodetic coordinates (latitude and longitude) to
     * Polar Stereographic coordinates (easting and northing), according to the current ellipsoid and Polar
     * Stereographic projection parameters. If any errors occur, error code(s) are returned by the function, otherwise
     * POLAR_NO_ERROR is returned.
     *
     * @param latitude  latitude, in radians
     * @param longitude longitude, in radians
     * @return error code
     * @see <a href="https://github.com/NASAWorldWind/WorldWindJava/blob/develop/src/gov/nasa/worldwind/geom/coords/PolarCoordConverter.java">gov.nasa.worldwind.geom.coords</a>
     */
    @SuppressWarnings("squid:S3776")
    public long convertGeodeticToPolarStereographic(double latitude, double longitude) {
        double dlam;
        double slat;
        double essin;
        double t;
        double rho;
        double powEs;
        long errorCode = POLAR_NO_ERROR;

        if ((latitude < -PI_OVER_2) || (latitude > PI_OVER_2)) {   /* latitude out of range */
            errorCode |= POLAR_LAT_ERROR;
        }
        if ((latitude < 0) && (southernHemisphere == 0)) {   /* latitude and Origin latitude in different hemispheres */
            errorCode |= POLAR_LAT_ERROR;
        }
        if ((latitude > 0) && (southernHemisphere == 1)) {   /* latitude and Origin latitude in different hemispheres */
            errorCode |= POLAR_LAT_ERROR;
        }
        if ((longitude < -PI) || (longitude > TWO_PI)) {  /* longitude out of range */
            errorCode |= POLAR_LON_ERROR;
        }

        if (errorCode == POLAR_NO_ERROR) {  /* no errors */

            if (Math.abs(Math.abs(latitude) - PI_OVER_2) < 1.0e-10) {
                easting = 0.0;
                northing = 0.0;
            } else {
                if (southernHemisphere != 0) {
                    longitude *= -1.0;
                    latitude *= -1.0;
                }
                dlam = longitude - polarOriginLong;
                if (dlam > PI) {
                    dlam -= TWO_PI;
                }
                if (dlam < -PI) {
                    dlam += TWO_PI;
                }
                slat = Math.sin(latitude);
                essin = es * slat;
                powEs = Math.pow((1.0 - essin) / (1.0 + essin), esOver2);
                t = Math.tan(PI_OVER_4 - latitude / 2.0) / powEs;

                if (Math.abs(Math.abs(polarOriginLat) - PI_OVER_2) > 1.0e-10)
                    rho = polarAMc * t / tc;
                else
                    rho = twoPolarA * t / e4;

                if (southernHemisphere != 0) { // TODO: What?
                    easting = -(rho * Math.sin(dlam) - polarFalseEasting);
                    northing = rho * Math.cos(dlam) + polarFalseNorthing;
                } else
                    easting = rho * Math.sin(dlam) + polarFalseEasting;
                northing = -rho * Math.cos(dlam) + polarFalseNorthing;
            }
        }
        return (errorCode);
    }

    public double getEasting() {
        return easting;
    }

    public double getNorthing() {
        return northing;
    }

    /**
     * The function Convert_Polar_Stereographic_To_Geodetic converts Polar
     * Stereographic coordinates (easting and northing) to geodetic
     * coordinates (latitude and longitude) according to the current ellipsoid
     * and Polar Stereographic projection Parameters. If any errors occur, the
     * code(s) are returned by the function, otherwise POLAR_NO_ERROR
     * is returned.
     *
     * @param easting  easting (X), in meters
     * @param northing northing (Y), in meters
     * @return error code
     */
    @SuppressWarnings("squid:S3776")
    public long convertPolarStereographicToGeodetic(double easting, double northing) {
        double dy = 0;
        double dx = 0;
        double rho = 0;
        double t;
        double phi;
        double sinPhi;
        double tempPHI = 0.0;
        double essin;
        double powEs;
        double deltaRadius;
        long errorCode = POLAR_NO_ERROR;
        double minEasting = polarFalseEasting - polarDeltaEasting;
        double maxEasting = polarFalseEasting + polarDeltaEasting;
        double minNorthing = polarFalseNorthing - polarDeltaNorthing;
        double maxNorthing = polarFalseNorthing + polarDeltaNorthing;

        if (easting > maxEasting || easting < minEasting) { /* easting out of range */
            errorCode |= POLAR_EASTING_ERROR;
        }
        if (northing > maxNorthing || northing < minNorthing) { /* northing out of range */
            errorCode |= POLAR_NORTHING_ERROR;
        }

        if (errorCode == POLAR_NO_ERROR) {
            dy = northing - polarFalseNorthing;
            dx = easting - polarFalseEasting;

            /* Radius of point with origin of false easting, false northing */
            rho = Math.sqrt(dx * dx + dy * dy);

            deltaRadius = Math.sqrt(polarDeltaEasting * polarDeltaEasting + polarDeltaNorthing * polarDeltaNorthing);

            if (rho > deltaRadius) { /* Point is outside of projection area */
                errorCode |= POLAR_RADIUS_ERROR;
            }
        }

        if (errorCode == POLAR_NO_ERROR) { /* no errors */
            if ((dy == 0.0) && (dx == 0.0)) {
                latitude = PI_OVER_2;
                longitude = polarOriginLong;

            } else {
                if (southernHemisphere != 0) {
                    dy *= -1.0;
                    dx *= -1.0;
                }

                if (Math.abs(Math.abs(polarOriginLat) - PI_OVER_2) > 1.0e-10)
                    t = rho * tc / (polarAMc);
                else
                    t = rho * e4 / (twoPolarA);
                phi = PI_OVER_2 - 2.0 * Math.atan(t);
                while (Math.abs(phi - tempPHI) > 1.0e-10) {
                    tempPHI = phi;
                    sinPhi = Math.sin(phi);
                    essin = es * sinPhi;
                    powEs = Math.pow((1.0 - essin) / (1.0 + essin), esOver2);
                    phi = PI_OVER_2 - 2.0 * Math.atan(t * powEs);
                }
                latitude = phi;
                longitude = polarOriginLong + Math.atan2(dx, -dy);

                if (longitude > PI)
                    longitude -= TWO_PI;
                else if (longitude < -PI)
                    longitude += TWO_PI;


                if (latitude > PI_OVER_2)  /* force distorted values to 90, -90 degrees */
                    latitude = PI_OVER_2;
                else if (latitude < -PI_OVER_2)
                    latitude = -PI_OVER_2;

                if (longitude > PI)  /* force distorted values to 180, -180 degrees */
                    longitude = PI;
                else if (longitude < -PI)
                    longitude = -PI;

            }
            if (southernHemisphere != 0) {
                latitude *= -1.0;
                longitude *= -1.0;
            }

        }
        return (errorCode);
    }

    /**
     * @return Latitude in radians.
     */
    public double getLatitude() {
        return latitude;
    }

    /**
     * @return Longitude in radians.
     */
    public double getLongitude() {
        return longitude;
    }

    /**
     * Computes the geodetic latitude given a spherical latitude in radians and a height h
     *
     * @param gcLat the geocentric latitude in radians
     * @param h     the height in meters
     * @return double the geodetic latitude
     **/
    public double gcLat2gdLat(double gcLat, double h) {

        double a = Math.pow(Math.sin(gcLat), 2);
        double rn = EARTH_RADIUS_EQ_M / (Math.sqrt(1 - WGS84_E2 * a));
        double b = 1 - WGS84_E2 * (rn / (rn + h));
        return Math.atan(Math.tan(gcLat) / b);

    }

    /**
     * Computes the geodetic latitude given a spherical latitude for a point in the surface of a sphere
     *
     * @param gcLat the geocentric latitude in radians
     * @return double the geodetic latitude
     **/
    public double gcLat2gdLat(double gcLat) {
        return Math.atan2(Math.tan(gcLat), (1 - WGS84_E2));
    }

    /**
     * Computes the geodetic latitude given a spherical latitude for a point in the surface of a sphere
     *
     * @param gcLat the geocentric latitude in degrees
     * @return double the geodetic latitude
     **/
    public double geodeticLatitudeToSphericalLatitudeDegrees(double gcLat) { // TODO normalize either degrees or radians usage
        double gcLatRad = Math.toRadians(gcLat);
        return Math.toDegrees(gcLat2gdLat(gcLatRad));
    }

    /**
     * Computes the geocentric latitude given a geodetic latitude for a point in the surface of an ellipsoid
     *
     * @param lat the geodetic latitude in radians
     * @return double the geocentric latitude in radians
     **/
    public double geodeticLatitudeToSphericalLatitudeRadians(double lat) { // TODO normalize either degrees or radians usage
        return Math.atan((1 - WGS84_E2) * Math.tan(lat));
    }

    /**
     * Computes the geocentric latitude given a geodetic latitude for a point in the surface of an ellipsoid
     *
     * @param gdLat the geodetic latitude in degrees
     * @return double the geocentric latitude
     **/
    public double gdLat2gcLatD(double gdLat) { // TODO normalize either degrees or radians usage
        double gdLatRad = Math.toRadians(gdLat);
        return Math.toDegrees(geodeticLatitudeToSphericalLatitudeRadians(gdLatRad));
    }

    public double conformal2latD(double confLat) {

        return Math.toDegrees(conformal2lat(Math.toRadians(confLat)));

    }

    public double conformal2lat(double confLat) {

        return confLat + (Math.pow(WGS84_E, 2) / 2 + 5 * Math.pow(WGS84_E, 4) / 24 + 1 * Math.pow(WGS84_E, 6) / 12 + 13 * Math.pow(WGS84_E, 8) / 360) * Math.sin(2 * confLat)
                + (7 * Math.pow(WGS84_E, 4) / 48 + 29 * Math.pow(WGS84_E, 6) / 240 + 811 * Math.pow(WGS84_E, 8) / 11520) * Math.sin(4 * confLat)
                + (7 * Math.pow(WGS84_E, 6) / 120 + 81 * Math.pow(WGS84_E, 8) / 1120) * Math.sin(6 * confLat)
                + (4279 * Math.pow(WGS84_E, 8) / 161280) * Math.sin(8 * confLat);

    }

    public double lat2ConformalD(double confLat) {

        return Math.toDegrees(lat2conformal(Math.toRadians(confLat)));

    }

    public double lat2conformal(double sphericalLat) {

        return sphericalLat + (Math.pow(WGS84_E, 2) / 2 + 5 * Math.pow(WGS84_E, 4) / 24 + 3 * Math.pow(WGS84_E, 6) / 32 + 281 * Math.pow(WGS84_E, 8) / 5760) * Math.sin(2 * sphericalLat)
                + (5 * Math.pow(WGS84_E, 4) / 48 + 7 * Math.pow(WGS84_E, 6) / 80 + 697 * Math.pow(WGS84_E, 8) / 11520) * Math.sin(4 * sphericalLat)
                + (13 * Math.pow(WGS84_E, 6) / 480 + 461 * Math.pow(WGS84_E, 8) / 13440) * Math.sin(6 * sphericalLat)
                + (1237 * Math.pow(WGS84_E, 8) / 161280) * Math.sin(8 * sphericalLat);

    }

}
