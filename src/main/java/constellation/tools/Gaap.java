package constellation.tools;

import com.menecats.polybool.Epsilon;
import com.menecats.polybool.PolyBool;
import com.menecats.polybool.models.Polygon;
import constellation.tools.geometry.AAP;
import constellation.tools.geometry.ConstellationSSPs;
import constellation.tools.geometry.FOV;
import constellation.tools.geometry.Geo;
import constellation.tools.math.Combination;
import constellation.tools.math.Pair;
import constellation.tools.math.Transformations;
import constellation.tools.output.ReportGenerator;
import org.orekit.data.DataContext;
import org.orekit.data.DataProvidersManager;
import org.orekit.data.DirectoryCrawler;
import org.orekit.time.AbsoluteDate;
import satellite.tools.Simulation;
import satellite.tools.assets.entities.Satellite;
import satellite.tools.structures.Ephemeris;
import satellite.tools.utils.Log;
import satellite.tools.utils.Utils;

import java.awt.geom.Path2D;
import java.io.File;
import java.util.*;
import java.util.stream.Collectors;

import static com.menecats.polybool.helpers.PolyBoolHelper.epsilon;

/**
 * Geodetic Access Areas Propagator
 **/
public class Gaap {

    private final Properties properties = Utils.loadProperties("gaap.properties");

    private final String START_DATE = (String) properties.get("start_date");
    private final String END_DATE = (String) properties.get("end_date");
    private final double TIME_STEP = Double.parseDouble((String) properties.get("time_step"));
    private final String OUTPUT_PATH = (String) properties.get("output_path");
    private final String SATELLITES_FILE = (String) properties.get("satellites_file");
    private final boolean DEBUG = Boolean.parseBoolean((String) properties.get("debug_mode"));
    private final double VISIBILITY_THRESHOLD = Double.parseDouble((String) properties.get("visibility_threshold"));
    private final double POLYGON_SEGMENTS = Double.parseDouble((String) properties.get("polygon_segments"));
    private final int MAX_SUBSET_SIZE = Integer.parseInt((String) properties.get("max_subset_size"));
    private final boolean USE_CONFORMAL_LATITUDE = Boolean.parseBoolean((String) properties.get("use_conformal_latitude"));

    private final List<Satellite> satelliteList = Utils.satellitesFromFile(SATELLITES_FILE);
    private final List<String> statistics = new ArrayList<>();
    private final List<ConstellationSSPs> constellationSSPs = new ArrayList<>();

    private Map<Double, Double> radii = new HashMap<>();

    ReportGenerator reportGenerator = new ReportGenerator(OUTPUT_PATH);

    public static void main(String[] args) {

        Gaap gaap = new Gaap();
        gaap.run();

    }

    public Gaap() {



    }

    public void run() {

        // configure Orekit data context
        var orekitData = new File("src/main/resources/orekit-data");
        if (!orekitData.exists()) {
            Log.fatal("Failed to find orekit-data folder " + orekitData.getAbsolutePath());
            System.exit(1);
        }
        DataProvidersManager manager = DataContext.getDefault().getDataProvidersManager();
        manager.addProvider(new DirectoryCrawler(orekitData));

        AbsoluteDate endDate = Utils.stamp2AD(END_DATE);
        AbsoluteDate startDate = Utils.stamp2AD(START_DATE);
        AbsoluteDate pointerDate = startDate;
        double scenarioDuration = endDate.durationFrom(startDate);

        // We compute a Utility "List of Lists", containing all possible overlapping combinations between regions.
        Combination comb = new Combination(satelliteList.size(), MAX_SUBSET_SIZE);
        final List<List<Integer>> combinationsList = comb.computeCombinations();

        List<AAP> nonEuclideanAAPs = new ArrayList<>();
        List<AAP> euclideanAAPs = new ArrayList<>();

        double lambdaMax = getLambdaMax(satelliteList.get(0).getElement("a"), VISIBILITY_THRESHOLD); // FIXME do I use this?

        while (pointerDate.compareTo(endDate) <= 0) {

            long timeSinceStart = Utils.stamp2unix(pointerDate.toString()) - Utils.stamp2unix(START_DATE);
            updateProgressBar(pointerDate.durationFrom(startDate), scenarioDuration);

            // Obtain the starting non-euclidean FOVs and their surface value
            List<FOV> nonEuclideanFOVs = computeStartingFOVsAt(satelliteList, pointerDate);

            // Store the SSPs (per Guido's request)
            List<Pair> SSPs = new ArrayList<>();
            nonEuclideanFOVs.forEach(FOV -> SSPs.add(new Pair(FOV.getReferenceLat(), FOV.getReferenceLon())));
            constellationSSPs.add(new ConstellationSSPs(timeSinceStart, SSPs));

            // Accumulated areas by number of satellites in visibility is stored in this array (idx = number of sats, value = area) // FIXME remove once 1.1 is implemented
            Map<Integer, Double> accumulatedAreas = new HashMap<>(MAX_SUBSET_SIZE);

            double surfaceInKm = 0D;

            for (List<Integer> combination : combinationsList) {

                List<Pair> nonEuclideanCoordinates = new ArrayList<>();
                List<Pair> euclideanCoordinates = new ArrayList<>();

                List<FOV> FOVsToIntersect = new ArrayList<>();
                combination.forEach(regionIndex -> FOVsToIntersect.add(nonEuclideanFOVs.get(regionIndex)));
                int poleProximity = checkPoleInclusion(FOVsToIntersect, lambdaMax);

                double referenceLat = poleProximity * 90; // FOVsToIntersect.get(0).getReferenceLat();
                double referenceLon = 0; // FOVsToIntersect.get(0).getReferenceLon();

                // If this is the starting FOV
                if (combination.size() <= 1) {
                    int fovIdx = combination.get(0);
                    FOV neFov = nonEuclideanFOVs.get(fovIdx);
                    surfaceInKm = neFov.getSurfaceKm2();
                    nonEuclideanCoordinates = Transformations.polygon2pairList(nonEuclideanFOVs.get(fovIdx).getPolygon());
                    euclideanCoordinates = Transformations.polygon2pairList(Transformations.toEuclideanPlane(neFov.getPolygon(), referenceLat, referenceLon));

                } else if (checkDistances(combination, nonEuclideanFOVs, lambdaMax)) {    // FIXME optimize transforming every starting AAP to stereographic

//                  combination.forEach(regionIndex -> FOVsToIntersect.add(euclideanFOVs.get(regionIndex)));

                    List<Path2D.Double> polygonsToIntersect = new ArrayList<>();
                    FOVsToIntersect.forEach(FOV -> polygonsToIntersect.add(Transformations.toEuclideanPlane(FOV.getPolygon(), referenceLat, referenceLon)));

                    Path2D.Double resultingPolygon =  new Path2D.Double(polygonsToIntersect.get(0));

                    for (Path2D.Double polygon : polygonsToIntersect) {
                        if (polygonsToIntersect.indexOf(polygon) == 0) continue;

                        resultingPolygon = intersectAndGetPolygon(resultingPolygon, polygon);

                    }

                    euclideanCoordinates = Transformations.polygon2pairList(resultingPolygon);

                    // Fixme use just polygon, etymologically "area" brings too much confusion

                    // When going back to the non euclidean plane we cannot use area2pairList since, naturally, the Area object cannot be filled in that plane
                    Path2D.Double nonEuclideanIntersection = Transformations.toNonEuclideanPlane(resultingPolygon, referenceLat, referenceLon);
                    nonEuclideanCoordinates = Transformations.polygon2pairList(nonEuclideanIntersection);

                    surfaceInKm = Geo.computeNonEuclideanSurface(nonEuclideanCoordinates) * 1E-6;

                } else {
                    surfaceInKm = 0D;
//                    if (SAVE_EMPTY_AREAS) accessAreas.put(combination.toString(), 0D);
                }

                // Area accumulator and store // 1.1 FIXME Replace with post-surface filter and accumulator
                if (accumulatedAreas.containsKey(combination.size())) {
                    accumulatedAreas.put(combination.size(), accumulatedAreas.get(combination.size()) + surfaceInKm);
                } else {
                    accumulatedAreas.put(combination.size(), surfaceInKm);
                }

                // Save AAPs
                nonEuclideanAAPs.add(new AAP(timeSinceStart, combination.size(), combination,
                        nonEuclideanCoordinates.stream().map(pair -> pair.lat).collect(Collectors.toList()),
                        nonEuclideanCoordinates.stream().map(pair -> pair.lon).collect(Collectors.toList()),
                        surfaceInKm));

                euclideanAAPs.add(new AAP(timeSinceStart, combination.size(), combination,
                        euclideanCoordinates.stream().map(pair -> pair.lat).collect(Collectors.toList()),
                        euclideanCoordinates.stream().map(pair -> pair.lon).collect(Collectors.toList()),
                        surfaceInKm));

            }

            // Save the results // FIXME Define a result object
            statistics.add(stringifyResults(pointerDate, accumulatedAreas));

            // Advance to the next time step
            pointerDate = pointerDate.shiftedBy(TIME_STEP);

        }

        saveAccessRegions(nonEuclideanAAPs); // FIXME DEBUG
        saveAccessRegions2(euclideanAAPs);
        saveSSPs(constellationSSPs);

        reportGenerator.saveAsCSV(statistics, "stats");

    }

    private void saveAccessRegions(List<AAP> AAPs) {

        for (int nOfGw = 1; nOfGw <= MAX_SUBSET_SIZE; nOfGw++) {
            int finalNOfGw = nOfGw;

            reportGenerator.saveAsJSON(AAPs.stream()
                            .filter(AAP -> AAP.getnOfGwsInSight() == finalNOfGw)
                            .collect(Collectors.toList()),"NEPolygons_" + nOfGw);
        }

        // DEBUG
        reportGenerator.saveAsJSON(AAPs.stream()
                .filter(AAP -> AAP.getDate() == 780000)
                .collect(Collectors.toList()),"NEPolygonsDebug");

    }

    private void saveAccessRegions2(List<AAP> AAPs) {

        for (int nOfGw = 1; nOfGw <= MAX_SUBSET_SIZE; nOfGw++) {
            int finalNOfGw = nOfGw;
            reportGenerator.saveAsJSON(AAPs.stream()
                    .filter(AAP -> AAP.getnOfGwsInSight() == finalNOfGw)
                    .collect(Collectors.toList()),"EPolygons_" + nOfGw);
        }

        // DEBUG
        reportGenerator.saveAsJSON(AAPs.stream()
                .filter(AAP -> AAP.getDate() == 780000)
                .collect(Collectors.toList()),"EPolygonsDebug");

    }

    private void saveSSPs(List<ConstellationSSPs> constellationSSPs) {
        reportGenerator.saveAsJSON(constellationSSPs, "ConstellationSSPs");
    }

    private void updateProgressBar(double current, double total) {

        double progress = Math.round(current * 100 * 100.00 / total) / 100.00;

        // Progress bar
        for (int i = 0; i < (int) progress; i++) {
            System.out.print("\b");
        }
        for (int i = 1; i < (int) progress; i++) {
            System.out.print(":");
        }
        System.out.print(" " + (int) progress + " %");

    }

    // FIXME REPLACE WITH POST-ANALYSIS
    private String stringifyResults(AbsoluteDate pointerDate, Map<Integer, Double> accumulatedAreas) {

        StringBuilder sb = new StringBuilder();
        sb.append(Utils.stamp2unix(pointerDate.toString()) - Utils.stamp2unix(START_DATE));

        for (Integer key : accumulatedAreas.keySet()) {
            sb.append(",");
            sb.append(accumulatedAreas.get(key));
        }

        return sb.toString();

    }

    /**
     * This method takes the satellite list, propagates orbits to the specified date and computes the corresponding
     * access area or FOV polygon as a Region object for each one.
     *
     * @param satelliteList a list of Satellite objects
     * @param date          an AbsoluteDate object
     * @return a List of Regions
     **/
    public List<FOV> computeStartingFOVsAt(List<Satellite> satelliteList, AbsoluteDate date) {

        Simulation simulation = new Simulation();
        List<FOV> FOVList = new ArrayList<>();

        // Assume all satellites at the same height and take the first to get the needed size elements // FIXME In process, do not assume every satellite at the same height

        for (Satellite satellite : satelliteList) {
            simulation.setSatellite(satellite);
            Ephemeris ephemeris = simulation.computeSSPAndGetEphemeris(date);

            double lambdaMax = getLambdaMax(satellite.getElement("a"), VISIBILITY_THRESHOLD); // NOTE MOVED LAMBDA OVER HERE
            Path2D.Double polygon = drawAAP(lambdaMax, ephemeris.getLatitude(),
                    ephemeris.getLongitude(), POLYGON_SEGMENTS);

            FOV FOV = new FOV(satellite.getId(), ephemeris.getLatitude(), ephemeris.getLongitude(), polygon);
            FOV.setSurface(Geo.computeNonEuclideanSurface(polygon));
            FOVList.add(FOV);
        }

        return FOVList;

    }

    /**
     * This method checks the distance between each and all pairs of coordinates corresponding to the satellites SSPs,
     * and returns whether the intersection is empty or not
     *
     * @param assetsToCheck A List of Integer depicting the indexes of the satellites that are being checked
     * @param FOVList       A List of Region objects that will provide the SSP coordinates for the assets being checked
     * @param lambdaMax     The maximum Earth Central Angle
     * @return boolean whether the intersection is empty or not
     **/
    private boolean checkDistances(List<Integer> assetsToCheck, List<FOV> FOVList, double lambdaMax) {

        // First we need to generate the combination list for the pair of assets that need checking
        Combination combination = new Combination();

        List<List<Integer>> pairsToCheck = combination.computeCombinations(assetsToCheck, 2);

        for (List<Integer> pairToCheck : pairsToCheck) {

            int r1Idx = pairToCheck.get(0);
            int r2Idx = pairToCheck.get(1);
            double distance = Geo.computeGeodesic(FOVList.get(r1Idx), FOVList.get(r2Idx));

            if (distance >= 2 * lambdaMax) {
                return false;
            }

        }

        return true;
    }

    /**
     * Checks whether some FOV in the provided List contains any of Earth's poles
     *
     * @param regionsToIntersect A List of Regions to check
     * @param lambdaMax          The Maximum Earth Central Angle of the regions to check
     * @return 0 if no FOV contains either the north or South Pole, 1 if some FOV contains the North Pole,
     * -1 if some FOV contains the South Pole
     **/
    private int checkPoleInclusion(List<FOV> regionsToIntersect, double lambdaMax) {

        // check proximity poles
        for (FOV FOV : regionsToIntersect) {
            int proximity = checkPoleInclusion(FOV, lambdaMax);
            if (proximity != 0) {
                return proximity;
            }
        }
        return 0;
    }

    /**
     * Checks whether the provided FOV contains any of Earth's poles
     *
     * @param FOV       A List of Regions to check
     * @param lambdaMax The Maximum Earth Central Angle of the regions to check
     * @return 0 if the FOV does not contain either the north or South Pole, 1 if it contains the North Pole,
     * -1 if it contains the South Pole
     **/
    private int checkPoleInclusion(FOV FOV, double lambdaMax) {

        if (Geo.computeGeodesic(FOV.getReferenceLat(), FOV.getReferenceLon(), 90, 0) <= lambdaMax) {
            return 1;
        } else if (Geo.computeGeodesic(FOV.getReferenceLat(), FOV.getReferenceLon(), -90, 0) <= lambdaMax) {
            return -1;
        }
        return 0;
    }


    private Path2D.Double intersectAndGetPolygon(Path2D.Double polygon1, Path2D.Double polygon2) {

        /* Martinez-Rueda clipping algorithm - Java port by Menecats: https://github.com/Menecats/polybool-java
        * Paper: https://www.sciencedirect.com/science/article/pii/S0965997813000379
        * apparently, this one is significantly faster than other algorithms (maybe incl. JTS's one?)
        * and can better handle special cases (like self-intersecting polygons). It's also supposed to be "simple".
        * */

        // Transform polygon1 to Polygol Polygon
        List<List<double[]>> regions1 = new ArrayList<>();
        List<double[]> coordinates = new ArrayList<>();
        List<Pair> pairList = Transformations.polygon2pairList(polygon1);

        pairList.forEach(pair -> coordinates.add(pair.getPoint()));

        regions1.add(coordinates);
        Polygon polyA = new Polygon(regions1);

        List<List<double[]>> regions2 = new ArrayList<>();
        List<double[]> coordinates2 = new ArrayList<>();
        List<Pair> pairList2 = Transformations.polygon2pairList(polygon2);

        pairList2.forEach(pair -> coordinates2.add(pair.getPoint()));

        regions2.add(coordinates2);
        Polygon polyB = new Polygon(regions2);

        Epsilon eps = epsilon();

        Path2D.Double intersectionPolygon = new Path2D.Double();

        if (polyA.getRegions().get(0).size() >= 3 && polyB.getRegions().get(0).size() >= 3) {

            Polygon intersection = PolyBool.intersect(eps, polyA, polyB);

            if (!intersection.getRegions().isEmpty()) {
                intersection.getRegions().get(0).forEach(pair -> {
                    if (intersection.getRegions().get(0).indexOf(pair) != 0) {
                        intersectionPolygon.lineTo(pair[0], pair[1]);
                    } else {
                        intersectionPolygon.moveTo(pair[0], pair[1]);
                    }
                });
            }

        }

        return intersectionPolygon;

    }

    /**
     * This method takes the satellite list, propagates orbits to the specified date and computes each corresponding
     * access area or FOV polygon for each one
     *
     * @param satellite the satellite for which the polygon is to be calculated
     * @param date      a timestamp date String
     * @return List
     **/
    @Deprecated
    public Path2D.Double getStartingPolygonAt(Satellite satellite, String date) {

        Simulation simulation = new Simulation();
        AbsoluteDate absoluteDate = Utils.stamp2AD(date);

        // Assume all satellites at the same height and take the first to get the needed size elements
        double lambdaMax = getLambdaMax(satelliteList.get(0).getElement("a"), VISIBILITY_THRESHOLD);
        simulation.setSatellite(satellite);
        Ephemeris ephemeris = simulation.computeSSPAndGetEphemeris(absoluteDate);
        return drawAAP(lambdaMax, ephemeris.getLatitude(),
                ephemeris.getLongitude(), POLYGON_SEGMENTS);

    }

    /**
     * This method returns a polygon that approximates a circular access area centered at (centerLat, centerLon) with
     * a radius of lambdaMax and the passed number of segments
     *
     * @param lambdaMax the Maximum Earth Central Angle
     * @param centerLat the latitude of the circle's center
     * @param centerLon the longitude of the circle's center
     * @param segments  the amount of segments for the polygon
     * @return a Path2D.Double containing the polygon (counter clock-wise direction)
     **/
    public Path2D.Double drawAAP(double lambdaMax, double centerLat, double centerLon, double segments) {

        // java.awt.geom
        Path2D.Double euclideanPolygon = new Path2D.Double();

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

            if (segment == 1) {
                euclideanPolygon.moveTo(pointerLat, pointerLon);
            } else {
                euclideanPolygon.lineTo(pointerLat, pointerLon);
            }

        }

        return euclideanPolygon;

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


}