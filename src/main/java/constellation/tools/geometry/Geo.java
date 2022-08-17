package constellation.tools.geometry;

import constellation.tools.math.Pair;
import net.sf.geographiclib.*;
import org.apache.commons.math3.geometry.euclidean.threed.Vector3D;
import satellite.tools.utils.Log;
import satellite.tools.utils.Utils;

import java.awt.geom.Path2D;
import java.awt.geom.PathIterator;
import java.io.*;
import java.util.ArrayList;
import java.util.List;

/**
 * Contains Geographic operations
 */
public class Geo {

    public static final double EARTH_RADIUS_EQ_M = 6378135.0D;
    public static final double WGS84_EQ_RADIUS_M = 6378137.0D;
    public static final double WGS84_EQ_RADIUS_KM = 6378.137D;
    public final double WGS84_F = 1 / 298.257223563;
    public final double WGS84_E2 = 0.00669437999014;
    public final double WGS84_E = 0.081819190842613;

    private final long POLAR_NO_ERROR = 0x0000;
    private final long POLAR_LAT_ERROR = 0x0001;
    private final long POLAR_LON_ERROR = 0x0002;
    private final long POLAR_ORIGIN_LAT_ERROR = 0x0004;
    private final long POLAR_ORIGIN_LON_ERROR = 0x0008;
    public final long POLAR_EASTING_ERROR = 0x0010;
    public final long POLAR_NORTHING_ERROR = 0x0020;
    private final long POLAR_A_ERROR = 0x0040;
    private final long POLAR_INV_F_ERROR = 0x0080;
    public final long POLAR_RADIUS_ERROR = 0x0100;

    private final double PI = 3.14159265358979323;
    private final double PI_OVER_2 = PI / 2.0;
    private final double PI_Over_4 = PI / 4.0;
    private final double TWO_PI = 2.0 * PI;

    /* Ellipsoid Parameters, default to WGS 84  */
    private double Polar_a = 6378137.0;                    /* Semi-major axis of ellipsoid in meters  */
    private double Polar_f = 1 / 298.257223563;            /* Flattening of ellipsoid  */
    private double es = 0.08181919084262188000;            /* Eccentricity of ellipsoid    */
    private double es_OVER_2 = .040909595421311;           /* es / 2.0 */
    private double Southern_Hemisphere = 0;                /* Flag variable */
    private double mc = 1.0;
    private double tc = 1.0;
    private double e4 = 1.0033565552493;
    private double Polar_a_mc = 6378137.0;                 /* Polar_a * mc */
    private double two_Polar_a = 12756274.0;               /* 2.0 * Polar_a */

    /* Polar Stereographic projection Parameters */
    private double Polar_Origin_Lat = ((PI * 90) / 180);   /* Latitude of origin in radians */
    private double Polar_Origin_Long = 0.0;                /* Longitude of origin in radians */
    private double Polar_False_Easting = 0.0;              /* False easting in meters */
    private double Polar_False_Northing = 0.0;             /* False northing in meters */

    /* Maximum variance for easting and northing values for WGS 84. */
    private double Polar_Delta_Easting = 12713601.0;
    private double Polar_Delta_Northing = 12713601.0;

    private double Easting;
    private double Northing;
    private double Latitude;
    private double Longitude;

    public Geo() {

    }

    /**
     * Checks whether the provided FOV contains any of Earth's poles
     *
     * @param FOV       A List of Regions to check
     * @param lambdaMax The Maximum Earth Central Angle of the regions to check
     * @return 0 if the FOV does not contain either the north or South Pole, 1 if it contains the North Pole,
     * -1 if it contains the South Pole
     **/
    public int checkPoleInclusion(FOV FOV, double lambdaMax) {

        if (computeGeodesic(FOV.getSspLat(), FOV.getSspLon(), 90, 0) <= lambdaMax) {
            return 1;
        } else if (computeGeodesic(FOV.getSspLat(), FOV.getSspLon(), -90, 0) <= lambdaMax) {
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
    public double computeGeodesic(FOV r1, FOV r2) {
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
     * vertices. This method uses the net.sf.geographiclib library.
     *
     * @param pairList A List containing the coordinates of the polygon
     * @return Double the computed area in meters squared
     **/
    public double computeNonEuclideanSurface(List<Pair> pairList) { // TODO: REMOVE IF 2 WORKS

        PolygonArea polygonArea = new PolygonArea(Geodesic.WGS84, false);

        for (Pair pair : pairList) {
            polygonArea.AddPoint(pair.lat, pair.lon);
        }

        PolygonResult result = polygonArea.Compute();
        return Math.abs(result.area);

    }

    /**
     * This method computes the geodetic area of a list of coordinates, given by Pair objects depicting the polygon's
     * vertices. This method uses the net.sf.geographiclib library. Paper: https://doi.org/10.1007/s00190-012-0578-z
     * Section 6 , C. F. F. Karney, Algorithms for geodesics, J. Geodesy 87, 43–55 (2013).
     *
     * @param pairList A List containing the coordinates of the polygon
     * @return Double the computed area in meters squared
     **/
    public double computeNonEuclideanSurface2(List<double[]> pairList) {

        PolygonArea polygonArea = new PolygonArea(Geodesic.WGS84, false);

        for (double[] pair : pairList) {
            polygonArea.AddPoint(pair[0], pair[1]);
        }

        PolygonResult result = polygonArea.Compute();
        return Math.abs(result.area);

    }

    // TODO change getLambda name and add a method without the threshold
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
     * This method returns a polygon that approximates a circular access area centered at (centerLat, centerLon) with
     * a radius of lambdaMax and the passed number of segments
     *
     * @param lambdaMax the Maximum Earth Central Angle
     * @param centerLat the latitude of the circle's center
     * @param centerLon the longitude of the circle's center
     * @param segments  the amount of segments for the polygon
     * @return a List of double[] containing the polygon (counter clock-wise direction)
     **/
    public List<double[]> drawCircularAAP(double lambdaMax, double centerLat, double centerLon, double segments) {

        List<double[]> coordinates = new ArrayList<>();

        double lambdaMaxRads = Math.toRadians(lambdaMax);
        double theta;

        // Obtain center latitude' in radians
        double centerLat_ = (Math.PI / 2.0) - (Math.toRadians(gdLat2gcLatD(centerLat)));

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

            // FIXME use the normalize angle of Geo netgeographiclib
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

//            double conformalLat = Transformations.lat2conformal(pointerLat);
//            conformalLat = pointerLat;
            coordinates.add(new double[]{gcLat2gdLatD(pointerLat), pointerLon});

        }
        return coordinates;
    }

    /** // TODO ship this to Utils
     * Reads a file containing asset(s) parameter(s) and returns a list of objects accordingly
     *
     * @return List<Pair>
     */
    public List<double[]> file2DoubleList(String fileName) {

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

    public int checkPoleInclusion(double[] coordinate, List<double[]> polygon) {

        // First transform the polygon into a Path2D
        Path2D.Double path2D = new Path2D.Double();
        PathIterator iterator = path2D.getPathIterator(null);

        int segment = 0;

        for (double[] pair : polygon) {

            if (segment == 0) {
                path2D.moveTo(pair[0], pair[1]);
            } else {
                path2D.lineTo(pair[0], pair[1]);
            }

            segment++;

        }

        // TODO FINISH THIS METHOD

        return 0;

    }

    public double[] computeAntipode(double lat, double lon) {

        double[] antipode = new double[2];
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
            euclideanPolygon.add(new double[]{this.Easting, this.Northing});

        }

        return euclideanPolygon;

    }

    /**
     * Takes the a pair of geographic coordinates and transforms them into the euclidean plane
     **/
    public static double[] toStereo(double lat, double lon, double referenceLatRads, double referenceLonRads) {

        double localRadius = Utils.EARTH_RADIUS_AVG_KM;

        double k = (2 * localRadius) / (1 + Math.sin(referenceLatRads) * Math.sin(lat) +
                Math.cos(referenceLatRads) * Math.cos(lat) * Math.cos(lon - referenceLonRads));

        double xStereo = k * Math.cos(lat) * Math.sin(lon - referenceLonRads);
        double yStereo = k * (Math.cos(referenceLatRads) * Math.sin(lat) - Math.sin(referenceLatRads) * Math.cos(lat) * Math.cos(lon - referenceLonRads));

        return new double[]{xStereo, yStereo};

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

//
//
//            double rho = Math.sqrt(Math.pow(xStereo, 2.000) + Math.pow(yStereo, 2.000));
//
//            double localRadius = Utils.EARTH_RADIUS_AVG_KM;
//
//            double c = 2 * Math.atan2(rho, 2.0 * localRadius);
//            double lat = Math.asin(Math.cos(c) * Math.sin(referenceLatRads) + (yStereo * Math.sin(c) * Math.cos(referenceLatRads)) / rho);
//            double lon;
//
//            // For exactly the poles, avoid indeterminate points in the equations
//            if (referenceLat == 90.0) {
//                lon = referenceLonRads + Math.atan2(xStereo, (-yStereo));
//            } else if (referenceLat == -90.0) {
//                lon = referenceLonRads + Math.atan2(xStereo, yStereo);
//            } else {
//                lon = referenceLonRads + Math.atan2((xStereo * Math.sin(c)), (rho * Math.cos(referenceLatRads) * Math.cos(c)
//                        - yStereo * Math.sin(referenceLatRads) * Math.sin(c)));
//            }
//
//            // Go back to degrees
//            lat = Math.toDegrees(lat);
//            lon = Math.toDegrees(lon);
//
//            while (lon < -180D) lon += 360;
//            while (lon > 180D) lon -= 360;

            GCSPolygon.add(new double[]{Math.toDegrees(this.Latitude), Math.toDegrees(this.Longitude)});

        }

        return GCSPolygon;

    }

    /**
     * The function setPolarStereographicParameters receives the ellipsoid parameters and Polar Stereograpic projection
     * parameters as inputs, and sets the corresponding state variables.  If any errors occur, error code(s) are
     * returned by the function, otherwise POLAR_NO_ERROR is returned.
     *
     * @param a                        Semi-major axis of ellipsoid, in meters
     * @param f                        Flattening of ellipsoid
     * @param Latitude_of_True_Scale   Latitude of true scale, in radians
     * @param Longitude_Down_from_Pole Longitude down from pole, in radians
     * @param False_Easting            Easting (X) at center of projection, in meters
     * @param False_Northing           Northing (Y) at center of projection, in meters
     * @return error code
     */
    public long setPolarStereographicParameters(double a, double f, double Latitude_of_True_Scale,
                                                double Longitude_Down_from_Pole, double False_Easting, double False_Northing)

    {
        double es2;
        double slat, clat;
        double essin;
        double one_PLUS_es, one_MINUS_es;
        double pow_es;
        double inv_f = 1 / f;
        final double epsilon = 1.0e-2;
        long Error_Code = POLAR_NO_ERROR;

        if (a <= 0.0)
        { /* Semi-major axis must be greater than zero */
            Error_Code |= POLAR_A_ERROR;
        }
        if ((inv_f < 250) || (inv_f > 350))
        { /* Inverse flattening must be between 250 and 350 */
            Error_Code |= POLAR_INV_F_ERROR;
        }
        if ((Latitude_of_True_Scale < -PI_OVER_2) || (Latitude_of_True_Scale > PI_OVER_2))
        { /* Origin Latitude out of range */
            Error_Code |= POLAR_ORIGIN_LAT_ERROR;
        }
        if ((Longitude_Down_from_Pole < -PI) || (Longitude_Down_from_Pole > TWO_PI))
        { /* Origin Longitude out of range */
            Error_Code |= POLAR_ORIGIN_LON_ERROR;
        }

        if (Error_Code == POLAR_NO_ERROR)
        { /* no errors */

            Polar_a = a;
            two_Polar_a = 2.0 * Polar_a;
            Polar_f = f;

            if (Longitude_Down_from_Pole > PI)
                Longitude_Down_from_Pole -= TWO_PI;
            if (Latitude_of_True_Scale < 0)
            {
                Southern_Hemisphere = 1;
                Polar_Origin_Lat = -Latitude_of_True_Scale;
                Polar_Origin_Long = -Longitude_Down_from_Pole;
            } else
            {
                Southern_Hemisphere = 0;
                Polar_Origin_Lat = Latitude_of_True_Scale;
                Polar_Origin_Long = Longitude_Down_from_Pole;
            }
            Polar_False_Easting = False_Easting;
            Polar_False_Northing = False_Northing;

            es2 = 2 * Polar_f - Polar_f * Polar_f;
            es = Math.sqrt(es2);
            es_OVER_2 = es / 2.0;

            if (Math.abs(Math.abs(Polar_Origin_Lat) - PI_OVER_2) > 1.0e-10)
            {
                slat = Math.sin(Polar_Origin_Lat);
                essin = es * slat;
                pow_es = Math.pow((1.0 - essin) / (1.0 + essin), es_OVER_2);
                clat = Math.cos(Polar_Origin_Lat);
                mc = clat / Math.sqrt(1.0 - essin * essin);
                Polar_a_mc = Polar_a * mc;
                tc = Math.tan(PI_Over_4 - Polar_Origin_Lat / 2.0) / pow_es;
            } else
            {
                one_PLUS_es = 1.0 + es;
                one_MINUS_es = 1.0 - es;
                e4 = Math.sqrt(Math.pow(one_PLUS_es, one_PLUS_es) * Math.pow(one_MINUS_es, one_MINUS_es));
            }
        }

        /* Calculate Radius */
        convertGeodeticToPolarStereographic(0, Polar_Origin_Long);

        Polar_Delta_Northing = Northing * 2; // Increased range for accepted easting and northing values
        Polar_Delta_Northing = Math.abs(Polar_Delta_Northing) + epsilon;
        Polar_Delta_Easting = Polar_Delta_Northing;

        return (Error_Code);
    }

    /**
     * The function Convert_Geodetic_To_Polar_Stereographic converts geodetic coordinates (latitude and longitude) to
     * Polar Stereographic coordinates (easting and northing), according to the current ellipsoid and Polar
     * Stereographic projection parameters. If any errors occur, error code(s) are returned by the function, otherwise
     * POLAR_NO_ERROR is returned.
     *
     * @param Latitude  latitude, in radians
     * @param Longitude Longitude, in radians
     * @return error code
     */
    public long convertGeodeticToPolarStereographic(double Latitude, double Longitude)
    {
        double dlam;
        double slat;
        double essin;
        double t;
        double rho;
        double pow_es;
        long Error_Code = POLAR_NO_ERROR;

        if ((Latitude < -PI_OVER_2) || (Latitude > PI_OVER_2))
        {   /* Latitude out of range */
            Error_Code |= POLAR_LAT_ERROR;
        }
        if ((Latitude < 0) && (Southern_Hemisphere == 0))
        {   /* Latitude and Origin Latitude in different hemispheres */
            Error_Code |= POLAR_LAT_ERROR;
        }
        if ((Latitude > 0) && (Southern_Hemisphere == 1))
        {   /* Latitude and Origin Latitude in different hemispheres */
            Error_Code |= POLAR_LAT_ERROR;
        }
        if ((Longitude < -PI) || (Longitude > TWO_PI))
        {  /* Longitude out of range */
            Error_Code |= POLAR_LON_ERROR;
        }

        if (Error_Code == POLAR_NO_ERROR)
        {  /* no errors */

            if (Math.abs(Math.abs(Latitude) - PI_OVER_2) < 1.0e-10)
            {
                Easting = 0.0;
                Northing = 0.0;
            } else
            {
                if (Southern_Hemisphere != 0)
                {
                    Longitude *= -1.0;
                    Latitude *= -1.0;
                }
                dlam = Longitude - Polar_Origin_Long;
                if (dlam > PI)
                {
                    dlam -= TWO_PI;
                }
                if (dlam < -PI)
                {
                    dlam += TWO_PI;
                }
                slat = Math.sin(Latitude);
                essin = es * slat;
                pow_es = Math.pow((1.0 - essin) / (1.0 + essin), es_OVER_2);
                t = Math.tan(PI_Over_4 - Latitude / 2.0) / pow_es;

                if (Math.abs(Math.abs(Polar_Origin_Lat) - PI_OVER_2) > 1.0e-10)
                    rho = Polar_a_mc * t / tc;
                else
                    rho = two_Polar_a * t / e4;


                if (Southern_Hemisphere != 0)
                {
                    Easting = -(rho * Math.sin(dlam) - Polar_False_Easting);
                    //Easting *= -1.0;
                    Northing = rho * Math.cos(dlam) + Polar_False_Northing;
                } else
                    Easting = rho * Math.sin(dlam) + Polar_False_Easting;
                Northing = -rho * Math.cos(dlam) + Polar_False_Northing;
            }
        }
        return (Error_Code);
    }

    public double getEasting()
    {
        return Easting;
    }

    public double getNorthing()
    {
        return Northing;
    }

    /**
     *  The function Convert_Polar_Stereographic_To_Geodetic converts Polar
     *  Stereographic coordinates (easting and northing) to geodetic
     *  coordinates (latitude and longitude) according to the current ellipsoid
     *  and Polar Stereographic projection Parameters. If any errors occur, the
     *  code(s) are returned by the function, otherwise POLAR_NO_ERROR
     *  is returned.
     *
     *  @param Easting Easting (X), in meters
     *  @param Northing Northing (Y), in meters
     *  @return error code
     */
    public long convertPolarStereographicToGeodetic (double Easting, double Northing)
    {
        double dy = 0, dx = 0;
        double rho = 0;
        double t;
        double PHI, sin_PHI;
        double tempPHI = 0.0;
        double essin;
        double pow_es;
        double delta_radius;
        long Error_Code = POLAR_NO_ERROR;
        double min_easting = Polar_False_Easting - Polar_Delta_Easting;
        double max_easting = Polar_False_Easting + Polar_Delta_Easting;
        double min_northing = Polar_False_Northing - Polar_Delta_Northing;
        double max_northing = Polar_False_Northing + Polar_Delta_Northing;

        if (Easting > max_easting || Easting < min_easting)
        { /* Easting out of range */
            Error_Code |= POLAR_EASTING_ERROR;
        }
        if (Northing > max_northing || Northing < min_northing)
        { /* Northing out of range */
            Error_Code |= POLAR_NORTHING_ERROR;
        }

        if (Error_Code == POLAR_NO_ERROR)
        {
            dy = Northing - Polar_False_Northing;
            dx = Easting - Polar_False_Easting;

            /* Radius of point with origin of false easting, false northing */
            rho = Math.sqrt(dx * dx + dy * dy);

            delta_radius = Math.sqrt(Polar_Delta_Easting * Polar_Delta_Easting + Polar_Delta_Northing * Polar_Delta_Northing);

            if(rho > delta_radius)
            { /* Point is outside of projection area */
                Error_Code |= POLAR_RADIUS_ERROR;
            }
        }

        if (Error_Code == POLAR_NO_ERROR)
        { /* no errors */
            if ((dy == 0.0) && (dx == 0.0))
            {
                Latitude = PI_OVER_2;
                Longitude = Polar_Origin_Long;

            }
            else
            {
                if (Southern_Hemisphere != 0)
                {
                    dy *= -1.0;
                    dx *= -1.0;
                }

                if (Math.abs(Math.abs(Polar_Origin_Lat) - PI_OVER_2) > 1.0e-10)
                    t = rho * tc / (Polar_a_mc);
                else
                    t = rho * e4 / (two_Polar_a);
                PHI = PI_OVER_2 - 2.0 * Math.atan(t);
                while (Math.abs(PHI - tempPHI) > 1.0e-10)
                {
                    tempPHI = PHI;
                    sin_PHI = Math.sin(PHI);
                    essin =  es * sin_PHI;
                    pow_es = Math.pow((1.0 - essin) / (1.0 + essin), es_OVER_2);
                    PHI = PI_OVER_2 - 2.0 * Math.atan(t * pow_es);
                }
                Latitude = PHI;
                Longitude = Polar_Origin_Long + Math.atan2(dx, -dy);

                if (Longitude > PI)
                    Longitude -= TWO_PI;
                else if (Longitude < -PI)
                    Longitude += TWO_PI;


                if (Latitude > PI_OVER_2)  /* force distorted values to 90, -90 degrees */
                    Latitude = PI_OVER_2;
                else if (Latitude < -PI_OVER_2)
                    Latitude = -PI_OVER_2;

                if (Longitude > PI)  /* force distorted values to 180, -180 degrees */
                    Longitude = PI;
                else if (Longitude < -PI)
                    Longitude = -PI;

            }
            if (Southern_Hemisphere != 0)
            {
                Latitude *= -1.0;
                Longitude *= -1.0;
            }

        }
        return (Error_Code);
    }

    /**
     * @return Latitude in radians.
     */
    public double getLatitude()
    {
        return Latitude;
    }

    /**
     * @return Longitude in radians.
     */
    public double getLongitude()
    {
        return Longitude;
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
    public double gcLat2gdLatD(double gcLat) { // TODO normalize either degrees or radians usage
        double gcLatRad = Math.toRadians(gcLat);
        return Math.toDegrees(gcLat2gdLat(gcLatRad));
    }

    /**
     * Computes the geocentric latitude given a geodetic latitude for a point in the surface of an ellipsoid
     *
     * @param lat the geodetic latitude in radians
     * @return double the geocentric latitude in radians
     **/
    public double gdLat2gcLat(double lat) { // TODO normalize either degrees or radians usage
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
        return Math.toDegrees(gdLat2gcLat(gdLatRad));
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
