package constellation.tools;

import com.menecats.polybool.Epsilon;
import com.menecats.polybool.PolyBool;
import com.menecats.polybool.models.Polygon;
import constellation.tools.geometry.AAP;
import constellation.tools.geometry.FOV;
import constellation.tools.geometry.Geo;
import constellation.tools.math.Combination;
import constellation.tools.math.Transformations;
import constellation.tools.reports.ReportGenerator;
import me.tongfei.progressbar.ProgressBar;
import me.tongfei.progressbar.ProgressBarStyle;
import org.orekit.data.DataContext;
import org.orekit.time.AbsoluteDate;
import satellite.tools.Simulation;
import satellite.tools.assets.entities.Satellite;
import satellite.tools.structures.Ephemeris;
import satellite.tools.utils.Log;
import satellite.tools.utils.Utils;

import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.*;
import java.util.stream.Collectors;

import static com.menecats.polybool.helpers.PolyBoolHelper.epsilon;
import static com.menecats.polybool.helpers.PolyBoolHelper.polygon;

/**
 * Dynamic Constellation Coverage Computer (D3CO)
 **/
public class D3CO {

    private static final Properties prop = Utils.loadProperties("config.properties");
    private static final String orekitPath = (String) prop.get("orekit_data_path");
    private final String START_DATE = (String) prop.get("start_date");
    private final long UNIX_START_DATE = stamp2unix(START_DATE);
    private final String END_DATE = (String) prop.get("end_date");
    private final double TIME_STEP = Double.parseDouble((String) prop.get("time_step"));
    private final String OUTPUT_PATH = (String) prop.get("output_path");
    private final String SATELLITES_FILE = (String) prop.get("satellites_file");
    private final String ROI_PATH = (String) prop.get("roi_path");
    private final boolean SAVE_SNAPSHOT = !((String) prop.get("snapshot")).isBlank();
    private final long SNAPSHOT = SAVE_SNAPSHOT ? Long.parseLong((String) prop.get("snapshot")) : 0L;
    private final boolean DEBUG_MODE = Boolean.parseBoolean((String) prop.get("debug_mode"));
    private final boolean SAVE_EUCLIDEAN = Boolean.parseBoolean((String) prop.get("save_euclidean"));
    private final boolean SAVE_GEOGRAPHIC = Boolean.parseBoolean((String) prop.get("save_geographic"));
    private final double VISIBILITY_THRESHOLD = Double.parseDouble((String) prop.get("visibility_threshold"));
    private final double POLYGON_SEGMENTS = Double.parseDouble((String) prop.get("polygon_segments"));
    private final int MAX_SUBSET_SIZE = Integer.parseInt((String) prop.get("max_subset_size"));

    private final List<Satellite> satelliteList = Utils.satellitesFromFile(SATELLITES_FILE);
    private final List<String> statistics = new ArrayList<>();

    ReportGenerator reportGenerator = new ReportGenerator(OUTPUT_PATH);
//    ProgressBar pb = new ProgressBar("Running", 100, ProgressBarStyle.ASCII);

    /**
     * Default constructor
     **/
    public D3CO() {


    }

    // FIXME URGENT: AM/PM!!!!!!!!!!!!
    public void run() {

        Simulation simulation = new Simulation(orekitPath);

        AbsoluteDate endDate = Utils.stamp2AD(END_DATE, DataContext.getDefault().getTimeScales().getUTC());
        AbsoluteDate startDate = Utils.stamp2AD(START_DATE);
        double scenarioDuration = endDate.durationFrom(startDate);

        // We compute a Utility "List of Lists", containing all possible overlapping combinations between regions.
        Combination comb = new Combination(satelliteList.size(), MAX_SUBSET_SIZE);
        final List<List<Integer>> combinationsList = comb.computeCombinations();

        List<AAP> nonEuclideanAAPs = new ArrayList<>();
        List<AAP> euclideanAAPs = new ArrayList<>();

        double lambdaMax = Geo.getLambdaMax(satelliteList.get(0).getElement("a"), VISIBILITY_THRESHOLD); // FIXME do I use this?

        if (DEBUG_MODE) Log.debug("Computing AAPs");

        long t0 = System.currentTimeMillis();

        try (ProgressBar pb = new ProgressBar("Obtaining AAPs", 100)) {

            pb.maxHint((long) scenarioDuration);
            for (AbsoluteDate t = startDate; t.compareTo(endDate) <= 0; t = t.shiftedBy(TIME_STEP)) {

                long timeSinceStart = stamp2unix(t.toString()) - stamp2unix(START_DATE);
                pb.stepTo((long) (Math.round(timeSinceStart * 100 * 100.00 / scenarioDuration) / 100.00));

                // Obtain the starting non-euclidean FOVs and their surface value
                List<FOV> nonEuclideanFOVs = computeFOVsAt(satelliteList, simulation, t);

                for (List<Integer> combination : combinationsList) {

                    List<double[]> nonEuclideanCoordinates = new ArrayList<>();
                    List<double[]> euclideanCoordinates = new ArrayList<>();

                    // assemble a list of the FOVs to be intersected at this time step:
                    List<FOV> FOVsToIntersect = new ArrayList<>();
                    combination.forEach(regionIndex -> FOVsToIntersect.add(nonEuclideanFOVs.get(regionIndex)));

                    // Get reference point for the projection //FIXME use satellite lambda
                    int poleProximity = checkPoleInclusion(FOVsToIntersect, lambdaMax);
                    double referenceLat = poleProximity * 90; // FOVsToIntersect.get(0).getReferenceLat();
                    double referenceLon = 0; // FOVsToIntersect.get(0).getReferenceLon();

                    // DEBUG
//                    if (timeSinceStart == SNAPSHOT) {
//                        Log.debug(referenceLat + "," + referenceLon);
//                        Log.debug("L: "+Geo.getLambdaMax(satelliteList.get(0).getElement("a"), VISIBILITY_THRESHOLD));
//                    }

                    // If this is the immediate FOV for a single satellite
                    if (combination.size() <= 1) {
                        int fovIdx = combination.get(0);
                        FOV neFov = nonEuclideanFOVs.get(fovIdx);
                        // nonEuclideanCoordinates = Transformations.doubleList2pairList(nonEuclideanFOVs.get(fovIdx).getPolygonCoordinates());
                        nonEuclideanCoordinates = nonEuclideanFOVs.get(fovIdx).getPolygonCoordinates();
                        euclideanCoordinates = Transformations.toEuclideanPlane(neFov.getPolygonCoordinates(),
                                referenceLat, referenceLon);

                    } else if (checkDistances(combination, nonEuclideanFOVs, lambdaMax)) {

                        List<List<double[]>> polygonsToIntersect = new ArrayList<>();
                        FOVsToIntersect.forEach(FOV -> polygonsToIntersect.add(Transformations.toEuclideanPlane(FOV.getPolygonCoordinates(), referenceLat, referenceLon)));

                        List<double[]> intersectedPolygon = new ArrayList<>(polygonsToIntersect.get(0));

                        // Obtain access polygon
                        for (List<double[]> polygon : polygonsToIntersect) {
                            if (polygonsToIntersect.indexOf(polygon) == 0) continue;
                            intersectedPolygon = intersectAndGetPolygon(intersectedPolygon, polygon);
                        }
                        // Intersected polygon euclidean coordinates
                        euclideanCoordinates = intersectedPolygon;

                        nonEuclideanCoordinates = Transformations.toNonEuclideanPlane(intersectedPolygon,
                                referenceLat, referenceLon);

                    }

                    // Save AAPs
                    nonEuclideanAAPs.add(new AAP(timeSinceStart, combination.size(), combination, nonEuclideanCoordinates,
                            nonEuclideanCoordinates.stream().map(pair -> pair[0]).collect(Collectors.toList()),
                            nonEuclideanCoordinates.stream().map(pair -> pair[1]).collect(Collectors.toList())));

                    if (SAVE_EUCLIDEAN) {
                        euclideanAAPs.add(new AAP(timeSinceStart, combination.size(), combination, euclideanCoordinates,
                                euclideanCoordinates.stream().map(pair -> pair[0]).collect(Collectors.toList()),
                                euclideanCoordinates.stream().map(pair -> pair[1]).collect(Collectors.toList())));
                    }

                }
            }
        }

        try (ProgressBar pb = new ProgressBar("Saving AAPs", 100)) {
            if (SAVE_GEOGRAPHIC) saveAAPs(nonEuclideanAAPs, "ne_polygons");
            pb.stepTo(25);
            if (SAVE_EUCLIDEAN) saveAAPs(euclideanAAPs, "e_polygons");
            pb.stepTo(50);
            if (SAVE_GEOGRAPHIC && SAVE_SNAPSHOT) saveAAPsAt(nonEuclideanAAPs, "snapshot_ne_polygons", SNAPSHOT);
            pb.stepTo(75);
            if (SAVE_EUCLIDEAN && SAVE_SNAPSHOT) saveAAPsAt(euclideanAAPs, "snapshot_e_polygons", SNAPSHOT);
            pb.stepTo(100);
        }

//        analyzeSurfaceCoverage(nonEuclideanAAPs);
        analyzeROICoverage(nonEuclideanAAPs);

        Log.info("Time to compute: " + (System.currentTimeMillis() - t0));

    }

    public void analyzeROICoverage(List<AAP> AAPs) {

        if (DEBUG_MODE) Log.debug("Computing ROI coverage");
        statistics.clear();

        // Load ROI Data:
        List<double[]> nonEuclideanROI = Geo.file2DoubleList(ROI_PATH);
        double roiSurface = Geo.computeNonEuclideanSurface2(nonEuclideanROI);

        // TODO: Generalize for any ROI

        double referenceLat = -90;
        double referenceLon = 0;

        List<double[]> euclideanROI = Transformations.toEuclideanPlane(nonEuclideanROI, referenceLat, referenceLon);

        // Timekeeping
        AbsoluteDate startDate = Utils.stamp2AD(START_DATE);
        AbsoluteDate pointerDate = Utils.stamp2AD(START_DATE);
        AbsoluteDate endDate = Utils.stamp2AD(END_DATE);
        double scenarioDuration = endDate.durationFrom(startDate);

        List<AAP> roiIntersections = new ArrayList<>();
        List<AAP> roiUnions = new ArrayList<>();

        try (ProgressBar pb = new ProgressBar("ROI coverage", 100)) {

            pb.maxHint((long) scenarioDuration);
            for (AbsoluteDate t = startDate; t.compareTo(endDate) <= 0; t = t.shiftedBy(TIME_STEP)) {

                long timeSinceStart = stamp2unix(t.toString()) - stamp2unix(START_DATE);
                pb.stepTo((long) (Math.round(timeSinceStart * 100 * 100.00 / scenarioDuration) / 100.00));

                long timeElapsed = stamp2unix(pointerDate.toString()) - stamp2unix(START_DATE);

                if (DEBUG_MODE) Log.debug(" t = " + pointerDate + " - unix = " + timeElapsed);

                // Group regions by number of satellites on sight, for this particular time step
                Map<Integer, List<AAP>> byAssetsInSight = mapByNOfAssets(AAPs, timeElapsed);

                double[] surfaceValues = new double[satelliteList.size()];

                // Perform intersection of AAPs with the ROI and surface area values calculation
                // For each number of assets
                byAssetsInSight.forEach((k, aaps) -> {

                    List<List<double[]>> unionQueue = new ArrayList<>();

                    // For each AAP with this number of assets in sight
                    aaps.forEach(aap -> {

                        // Intersections with ROI
                        List<double[]> eIntersection = intersectAndGetPolygon(euclideanROI,
                                Transformations.toEuclideanPlane(aap.getGeoCoordinates(),
                                        referenceLat, referenceLon));

                        if (eIntersection.size() >= 3) {
                            // Collections.reverse(eIntersection);
                            unionQueue.add(eIntersection);
                        }
                        List<double[]> neIntersection = Transformations.toNonEuclideanPlane(eIntersection, referenceLat, referenceLon);

                        AAP intersectionAAP = new AAP(timeElapsed, k, aap.getGwsInSight(), neIntersection,
                                neIntersection.stream().map(pair -> pair[0]).collect(Collectors.toList()),
                                neIntersection.stream().map(pair -> pair[1]).collect(Collectors.toList()));
                        roiIntersections.add(intersectionAAP);

                    });

                    if (!unionQueue.isEmpty()) {

                        Polygon union = polyUnion(unionQueue);

                        union.getRegions().forEach(region -> {
                            List<double[]> neIntersection = Transformations.toNonEuclideanPlane(region, referenceLat, referenceLon);
                            surfaceValues[k - 1] = surfaceValues[k - 1] + Geo.computeNonEuclideanSurface2(neIntersection);
                            AAP unionAAP = new AAP(timeElapsed, k, null, neIntersection,
                                    neIntersection.stream().map(pair -> pair[0]).collect(Collectors.toList()),
                                    neIntersection.stream().map(pair -> pair[1]).collect(Collectors.toList()));
                            roiUnions.add(unionAAP);
                        });
                    }

                });

                StringBuilder sb = new StringBuilder(timeElapsed + "");

                for (double surface : surfaceValues) {
                    sb.append(",");
                    double percentage = Math.round(((surface / roiSurface) * 100.00000) * 100000d) / 100000d;
                    sb.append(percentage);
                }

                statistics.add(sb.toString());

                // Advance timestep
                pointerDate = pointerDate.shiftedBy(TIME_STEP);

            }
        }

        if (SAVE_GEOGRAPHIC && SAVE_SNAPSHOT) {
            saveAAPsAt(roiIntersections, "snapshot_aaps_intersection", SNAPSHOT);
            saveAAPsAt(roiUnions, "snapshot_aaps_union", SNAPSHOT);
        }

        reportGenerator.saveAsCSV(statistics, "coverage");

    }

    /**
     * Performs the union of a list of polygons
     **/
    private Polygon polyUnion(List<List<double[]>> unionQueue) {

        Polygon union = new Polygon();

        // Union of all ROI intersections
        try {
            Epsilon eps = epsilon(0.0001);
            Polygon result = polygon(unionQueue.get(0));
            PolyBool.Segments segments = PolyBool.segments(eps, result);
            for (int i = 1; i < unionQueue.size(); i++) {
                PolyBool.Segments seg2 = PolyBool.segments(eps, polygon(unionQueue.get(i)));
                PolyBool.Combined comb = PolyBool.combine(eps, segments, seg2);
                segments = PolyBool.selectUnion(comb);
            }

            union = PolyBool.polygon(eps, segments);

        } catch (IndexOutOfBoundsException e1) {
            Log.error("IndexOutOfBoundsException " + e1.getMessage());
            Log.error("polygons to be or list size: " + unionQueue.size());
            unionQueue.forEach(region -> {
                Log.error("Region " + unionQueue.indexOf(region) + " size: " + unionQueue.size());
            });
        } catch (RuntimeException e2) {
            Log.error(e2.getMessage());
            Log.error("RuntimeException " + unionQueue.size());
            if (DEBUG_MODE) {
                unionQueue.forEach(region -> Log.error("Region " + unionQueue.indexOf(region) + " size: " + unionQueue.size()));
            }
        }

        return union;

    }

    private Map<Integer, List<AAP>> mapByNOfAssets(List<AAP> AAPs, long timeElapsed) {

        Map<Integer, List<AAP>> byAssetsInSight = new LinkedHashMap<>(MAX_SUBSET_SIZE);
        AAPs.stream().filter(AAP -> AAP.getDate() == timeElapsed).forEach(AAP -> {
            int nAssets = AAP.getnOfGwsInSight();
            // accumulatedAreas.putIfAbsent(nAssets, surfaceInKm2);
            if (byAssetsInSight.containsKey(nAssets)) {
                byAssetsInSight.get(nAssets).add(AAP);
            } else {
                List<AAP> aapList = new ArrayList<>();
                aapList.add(AAP);
                byAssetsInSight.put(nAssets, aapList);
            }
        });

        return byAssetsInSight;

    }

    public void analyzeSurfaceCoverage(List<AAP> AAPs) {

        statistics.clear();

        AbsoluteDate endDate = Utils.stamp2AD(END_DATE);
        AbsoluteDate pointerDate = Utils.stamp2AD(START_DATE);

        // Accumulated areas by number of satellites in visibility is stored in this array (idx = number of sats, value = area) // FIXME remove eventually
        Map<Integer, Double> accumulatedAreas = new HashMap<>(MAX_SUBSET_SIZE);

        while (pointerDate.compareTo(endDate) <= 0) {   // TODO while/filter can be improved

            accumulatedAreas.clear();

            long timeSinceStart = stamp2unix(pointerDate.toString()) - stamp2unix(START_DATE);
            AAPs.stream().filter(AAP -> AAP.getDate() == timeSinceStart).forEach(AAP -> {

                double surfaceInKm2 = Geo.computeNonEuclideanSurface2(AAP.getGeoCoordinates()) * 1E-6; // AAP.getSurfaceInKm2(); //
                int nAssets = AAP.getnOfGwsInSight();
                // accumulatedAreas.putIfAbsent(nAssets, surfaceInKm2);
                if (accumulatedAreas.containsKey(nAssets)) {
                    accumulatedAreas.put(nAssets, accumulatedAreas.get(nAssets) + surfaceInKm2);
                } else {
                    accumulatedAreas.put(nAssets, surfaceInKm2);
                }
                AAP.setSurfaceInKm2(surfaceInKm2);

            });

            // Save the results // FIXME Define a result object
            statistics.add(stringifyResults(pointerDate, accumulatedAreas));

            // Advance to the next time step
            pointerDate = pointerDate.shiftedBy(TIME_STEP);

        }

        reportGenerator.saveAsCSV(statistics, "stats");

    }

    private void saveAAPs(List<AAP> AAPs, String fileName) {

        for (int nOfGw = 1; nOfGw <= MAX_SUBSET_SIZE; nOfGw++) {
            int finalNOfGw = nOfGw;

            reportGenerator.saveAsJSON(AAPs.stream()
                    .filter(AAP -> AAP.getnOfGwsInSight() == finalNOfGw)
                    .collect(Collectors.toList()), fileName + "_" + nOfGw);
        }

    }

    private void saveAAPsAt(List<AAP> AAPs, String fileName, long time) {
        reportGenerator.saveAsJSON(AAPs.stream()
                .filter(AAP -> AAP.getDate() == time)
                .collect(Collectors.toList()), fileName);
    }

    // TODO REPLACE WITH POST-ANALYSIS
    private String stringifyResults(AbsoluteDate pointerDate, Map<Integer, Double> accumulatedAreas) {

        StringBuilder sb = new StringBuilder();
        sb.append(stamp2unix(pointerDate.toString()) - UNIX_START_DATE);

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
    private List<FOV> computeFOVsAt(List<Satellite> satelliteList, Simulation simulation, AbsoluteDate date) {

        List<FOV> FOVList = new ArrayList<>();

        for (Satellite satellite : satelliteList) {
            simulation.setSatellite(satellite);
            Ephemeris ephemeris = simulation.computeSSPAndGetEphemeris(date);


            double lambdaMax = Geo.getLambdaMax(satellite.getElement("a"), VISIBILITY_THRESHOLD);
            List<double[]> poly = Geo.drawCircularAAP(lambdaMax, ephemeris.getLatitude(), ephemeris.getLongitude(), POLYGON_SEGMENTS);

            FOV FOV = new FOV(satellite.getId(), ephemeris.getLatitude(), ephemeris.getLongitude(), poly);
            FOV.setPolygonCoordinates(poly);

            double surface = Geo.computeNonEuclideanSurface2(poly);

            FOV.setSurface(surface);
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

        for (List<Integer> pairToCheck : pairsToCheck) { // FIXME Use the biggest lambda

            int r1Idx = pairToCheck.get(0);
            int r2Idx = pairToCheck.get(1);

            // FIX? // TODO CHECK THIS
            double lambda1 = Geo.getLambdaMax(satelliteList.get(FOVList.get(r1Idx).getSatId()).getElement("a"), VISIBILITY_THRESHOLD);
            double lambda2 = Geo.getLambdaMax(satelliteList.get(FOVList.get(r2Idx).getSatId()).getElement("a"), VISIBILITY_THRESHOLD);
            double distance = Geo.computeGeodesic(FOVList.get(r1Idx), FOVList.get(r2Idx));

            if (distance >= (lambda1 + lambda2)) {
                return false;
            }

//            if (distance >= (2 * lambdaMax)) {
//                return false;
//            }

//            if (distance >= (2 * lambdaMax)) {
//                return false;
//            }

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

            // FIX // TODO REVISIT THIS
            double lambda = Geo.getLambdaMax(satelliteList.get(FOV.getSatId()).getElement("a"), VISIBILITY_THRESHOLD);

            int proximity = Geo.checkPoleInclusion(FOV, lambda);
            if (proximity != 0) {
                return proximity;
            }
        }
        return -1;

    }

    // TODO: this can be improved inheriting the library's intersection capabilities, for now, we dont trust them

    /**
     * This method takes two polygons, and returns their intersection using the Martinez-Rueda Algorithm.
     *
     * @see <a href="https://github.com/Menecats/polybool-java">Menecats-Polybool</a>
     * @see <a href="https://www.sciencedirect.com/science/article/pii/S0965997813000379">Martinez-Rueda clipping algorithm</a>
     **/
    private List<double[]> intersectAndGetPolygon(List<double[]> polygonA, List<double[]> polygonB) {

        List<List<double[]>> regions1 = new ArrayList<>();
        regions1.add(polygonA);

        List<List<double[]>> regions2 = new ArrayList<>();
        regions2.add(polygonB);

        Polygon polyA = new Polygon(regions1);
        Polygon polyB = new Polygon(regions2);
        Polygon intersection = new Polygon();

        if (polyA.getRegions().get(0).size() >= 3 && polyB.getRegions().get(0).size() >= 3) {
            Epsilon eps = epsilon();
            intersection = PolyBool.intersect(eps, polyA, polyB);
        }

        if (intersection.getRegions().size() > 0) {
            return intersection.getRegions().get(0);
        } else {
            return new ArrayList<>();
        }

    }

    // TODO: FIX URGENT IN SATELLITE TOOLS
    public static long stamp2unix(String timestamp) {
        SimpleDateFormat dateFormat = new SimpleDateFormat("yyyy-MM-dd'T'HH:mm:ss.SSS");
        dateFormat.setTimeZone(TimeZone.getTimeZone("UTCG"));
        Date parsedDate = new Date();

        try {
            parsedDate = dateFormat.parse(timestamp);
            return parsedDate.getTime();
        } catch (ParseException var4) {
            var4.printStackTrace();
            return parsedDate.getTime();
        }
    }


}