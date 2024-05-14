package constellation.tools;

import com.menecats.polybool.models.Polygon;
import constellation.tools.geometry.AccessAreaPolygon;
import constellation.tools.geometry.AccessRegion;
import constellation.tools.geometry.Geographic;
import constellation.tools.geometry.OblateAccessRegion;
import constellation.tools.math.Combination;
import constellation.tools.math.TimedMetricsRecord;
import constellation.tools.math.Transformations;
import constellation.tools.operations.PolygonOperator;
import constellation.tools.utilities.AppConfig;
import constellation.tools.utilities.Reports;
import constellation.tools.utilities.TimeUtils;
import me.tongfei.progressbar.ProgressBar;
import org.apache.commons.configuration2.ex.ConfigurationException;
import org.orekit.data.DataContext;
import org.orekit.frames.FramesFactory;
import org.orekit.time.AbsoluteDate;
import org.orekit.utils.IERSConventions;
import satellite.tools.Simulation;
import satellite.tools.assets.entities.Satellite;
import satellite.tools.structures.Ephemeris;
import satellite.tools.utils.Log;
import satellite.tools.utils.Utils;

import java.util.*;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.AtomicReference;
import java.util.stream.Collectors;

/**
 * Dynamic Constellation Coverage Computer (D3CO)
 **/
public class D3CO implements Runnable {

    private static final Properties prop = Utils.loadProperties("config.properties");
    private static final String orekitPath = (String) prop.get("orekit_data_path");
    private final String START_DATE = (String) prop.get("start_date");
    private final long UNIX_START_DATE = TimeUtils.stamp2unix(START_DATE);
    private final String END_DATE = (String) prop.get("end_date");
    private final double TIME_STEP = Double.parseDouble((String) prop.get("time_step"));
    private final String OUTPUT_PATH = (String) prop.get("output_path");
    private final String SATELLITES_FILE = (String) prop.get("satellites_file");
    private final String ROI_PATH = (String) prop.get("roi_path");
    private final String POSITIONS_PATH = (String) prop.get("positions_path");
    private final boolean SAVE_SNAPSHOT = !((String) prop.get("snapshot")).isBlank();
    private final long SNAPSHOT = SAVE_SNAPSHOT ? Long.parseLong((String) prop.get("snapshot")) : 0L;
    private final boolean DEBUG_MODE = Boolean.parseBoolean((String) prop.get("debug_mode"));
    private final boolean PROPAGATE_INTERNALLY = Boolean.parseBoolean((String) prop.get("propagate_internally"));
    private final boolean SAVE_EUCLIDEAN = Boolean.parseBoolean((String) prop.get("save_euclidean"));
    private final boolean SAVE_GEOGRAPHIC = Boolean.parseBoolean((String) prop.get("save_geographic"));
    private final double VISIBILITY_THRESHOLD = Double.parseDouble((String) prop.get("visibility_threshold"));
    private final int POLYGON_SEGMENTS = Integer.parseInt((String) prop.get("polygon_segments"));
    private final double POLYGON_EPSILON = Double.parseDouble((String) prop.get("polygon_epsilon"));
    private final double LAMBDA_EXCLUSION = Double.parseDouble((String) prop.get("lambda_exclusion"));
    private final int MAX_SUBSET_SIZE = Integer.parseInt((String) prop.get("max_subset_size"));

    private List<Satellite> satelliteList;
    private final List<String> statistics = new ArrayList<>();
    private List<Map<Long, Ephemeris>> constellation;

    Reports reports = new Reports(OUTPUT_PATH);

    private Geographic geographic = new Geographic();

    AppConfig appConfig;

    // TODO: Add ban list

    /**
     * Default constructor
     **/
    public D3CO() {
        try {
            appConfig = new AppConfig("config.properties");
        } catch (ConfigurationException e) {
            Log.warn("Error trying to load properties. ");
            System.exit(100);
        }
        satelliteList = Utils.satellitesFromFile(SATELLITES_FILE);

    }

    @Override
    public void run() {

        if (!PROPAGATE_INTERNALLY) {
            Log.debug("Reading positions from directory: " + POSITIONS_PATH);
            constellation = reports.positionsFromPath(POSITIONS_PATH, satelliteList.size());
        }

        Simulation simulation = new Simulation(orekitPath);
        simulation.setParams(START_DATE, END_DATE, TIME_STEP, VISIBILITY_THRESHOLD);
        simulation.setInertialFrame(FramesFactory.getEME2000());
        simulation.setFixedFrame(FramesFactory.getITRF(IERSConventions.IERS_2010, true));
        Log.info(simulation.getStartTime() + "-" + simulation.getEndTime());

        AbsoluteDate startDate = Utils.stamp2AD(START_DATE);
        AbsoluteDate endDate = Utils.stamp2AD(END_DATE, DataContext.getDefault().getTimeScales().getUTC());
        double scenarioDuration = endDate.durationFrom(startDate);

        // We compute a Utility "List of Lists", containing all possible overlapping combinations between regions.
        Combination comb = new Combination(satelliteList.size(), MAX_SUBSET_SIZE);
        final List<List<Integer>> combinationsList = comb.computeCombinations();

        Log.info("Satellites: " + '\n' + satelliteList.stream().map(s -> s.getElements().toString() + '\n').collect(Collectors.toList()) + '\n');

        List<AccessAreaPolygon> nonEuclideanAccessAreaPolygons = new ArrayList<>();
        List<AccessAreaPolygon> euclideanAccessAreaPolygons = new ArrayList<>();

        int avoided = 0;
        int performed = 0;
        int escaped = 0;

        try (ProgressBar pb = new ProgressBar("Obtaining AAPs", 100)) {

            pb.maxHint((long) scenarioDuration);

            for (AbsoluteDate t = startDate; t.compareTo(endDate) <= 0; t = t.shiftedBy(TIME_STEP)) {

                long timeElapsed = TimeUtils.stamp2unix(t.toString()) - TimeUtils.stamp2unix(START_DATE);
                pb.stepTo((long) (Math.round(timeElapsed * 100 * 100.00 / scenarioDuration) / 100.00));

                // Obtain the starting non-euclidean FOVs
                List<AccessRegion> nonEuclideanAccessRegions = computeAccessRegionsAt(satelliteList, simulation, t, timeElapsed);

                for (List<Integer> combination : combinationsList) {

                    List<double[]> nonEuclideanCoordinates = new ArrayList<>();
                    List<double[]> euclideanCoordinates = new ArrayList<>();

                    // assemble a list of the FOVs to be intersected at this time step:
                    List<AccessRegion> AccessRegionsToIntersect = new ArrayList<>();
                    combination.forEach(regionIndex -> AccessRegionsToIntersect.add(nonEuclideanAccessRegions.get(regionIndex)));

                    // Get reference point for the projection //FIXME use satellite lambda
                    double referenceLat = 0;
                    double referenceLon = 0;

                    // If this is the immediate FOV for a single satellite
                    if (combination.size() <= 1) {
                        int fovIdx = combination.get(0);
                        AccessRegion neAccessRegion = nonEuclideanAccessRegions.get(fovIdx);
                        nonEuclideanCoordinates = nonEuclideanAccessRegions.get(fovIdx).getPolygonCoordinates();
                        euclideanCoordinates = Transformations.toEuclideanPlane(neAccessRegion.getPolygonCoordinates(),
                                referenceLat, referenceLon);

                    } else if (checkDistances(combination, nonEuclideanAccessRegions) && !AccessRegionsToIntersect.isEmpty()) {

                        performed++;
                        List<List<double[]>> intersectionQueue = new LinkedList<>();
                        AccessRegionsToIntersect.forEach(AccessRegion -> intersectionQueue.add(Transformations.toEuclideanPlane(AccessRegion.getPolygonCoordinates(), referenceLat, referenceLon)));

                        // Obtain access polygon
                        PolygonOperator po = new PolygonOperator();
                        Polygon intersection = po.polyIntersect(intersectionQueue, POLYGON_EPSILON);
                        if (intersection.getRegions().size() > 0 && intersection.getRegions().get(0).size() > 2) {
                            if (SAVE_EUCLIDEAN) {
                                euclideanCoordinates = intersection.getRegions().get(0);
                            }
                            nonEuclideanCoordinates = Transformations.toNonEuclideanPlane(intersection.getRegions().get(0),
                                    referenceLat, referenceLon);
                        } else {
                            escaped++;
                        }

                    } else {
                        avoided++;
                    }

                    // Save AAPs
                    if (nonEuclideanCoordinates.size() > 2) {
                        nonEuclideanAccessAreaPolygons.add(new AccessAreaPolygon(timeElapsed, combination.size(), combination, nonEuclideanCoordinates,
                                nonEuclideanCoordinates.stream().map(pair -> pair[0]).collect(Collectors.toList()),
                                nonEuclideanCoordinates.stream().map(pair -> pair[1]).collect(Collectors.toList()),
                                referenceLat, referenceLon));
                    }

                    if (SAVE_EUCLIDEAN && euclideanCoordinates.size() > 0) {
                        euclideanAccessAreaPolygons.add(new AccessAreaPolygon(timeElapsed, combination.size(), combination, euclideanCoordinates,
                                euclideanCoordinates.stream().map(pair -> pair[0]).collect(Collectors.toList()),
                                euclideanCoordinates.stream().map(pair -> pair[1]).collect(Collectors.toList()),
                                referenceLat, referenceLon));
                    }

                }
            }
        }

        Log.info("Performed: " + performed + " - Avoided: " + avoided + " - Escaped: " + escaped);

        if (SAVE_GEOGRAPHIC || SAVE_EUCLIDEAN || SAVE_SNAPSHOT) {
            try (ProgressBar pb = new ProgressBar("Saving AAPs", 100)) {
                if (SAVE_GEOGRAPHIC) saveAAPs(nonEuclideanAccessAreaPolygons, "ne_polygons");
                pb.stepTo(25);
                if (SAVE_EUCLIDEAN) saveAAPs(euclideanAccessAreaPolygons, "e_polygons");
                pb.stepTo(50);
                if (SAVE_SNAPSHOT) saveAAPsAt(nonEuclideanAccessAreaPolygons, "snapshot_ne_polygons", SNAPSHOT);
                pb.stepTo(75);
                if (SAVE_SNAPSHOT && SAVE_EUCLIDEAN) saveAAPsAt(euclideanAccessAreaPolygons, "snapshot_e_polygons", SNAPSHOT);
                pb.stepTo(100);
            }
        }

        //        analyzeSurfaceCoverage(nonEuclideanAAPs);
        analyzeROICoverage(nonEuclideanAccessAreaPolygons);

    }

    public void analyzeROICoverage(List<AccessAreaPolygon> AccessAreaPolygons) {

        statistics.clear();

        // Load ROI Data:
        List<double[]> nonEuclideanROI = geographic.file2DoubleList(ROI_PATH);
        double roiSurface = geographic.computeNonEuclideanSurface(nonEuclideanROI);
        Log.info("ROI Surface: " + roiSurface);

        // Timekeeping
        AbsoluteDate startDate = Utils.stamp2AD(START_DATE);
        AbsoluteDate endDate = Utils.stamp2AD(END_DATE);
        double scenarioDuration = endDate.durationFrom(startDate);
        long startTimestamp = TimeUtils.stamp2unix(START_DATE);

        List<AccessAreaPolygon> roiIntersections = new ArrayList<>();
        List<AccessAreaPolygon> roiUnions = new ArrayList<>();
        List<TimedMetricsRecord> timeSeriesData = new ArrayList<>();
        AtomicInteger changes = new AtomicInteger();

        try (ProgressBar pb = new ProgressBar("ROI coverage", 100)) {

            pb.maxHint((long) scenarioDuration);
            for (AbsoluteDate t = startDate; t.compareTo(endDate) <= 0; t = t.shiftedBy(TIME_STEP)) {

                long timeElapsed = TimeUtils.stamp2unix(t.toString()) - startTimestamp;
                pb.stepTo((long) (Math.round(timeElapsed * 100 * 100.00 / scenarioDuration) / 100.00));

                // Group regions by number of satellites on sight, for this particular time step
                Map<Integer, List<AccessAreaPolygon>> byAssetsInSight = mapByNOfAssets(AccessAreaPolygons, timeElapsed);

                double[] surfaceValues = new double[satelliteList.size()];

                // Starting euclidean ROI
                AtomicReference<Double> roiReferenceLat = new AtomicReference<>(byAssetsInSight.get(1).get(0).getReferenceLat());
                AtomicReference<Double> roiReferenceLon = new AtomicReference<>(byAssetsInSight.get(1).get(0).getReferenceLon());
                AtomicReference<List<double[]>> euclideanROI = new AtomicReference<>(Transformations.toEuclideanPlane(nonEuclideanROI, roiReferenceLat.get(), roiReferenceLon.get()));

                TimedMetricsRecord timedMetricsRecord = new TimedMetricsRecord(timeElapsed, satelliteList.size());

                // Perform intersection of AAPs with the ROI and surface area values calculation
                // For each number of assets
                byAssetsInSight.forEach((k, aaps) -> {

                    // TODO: Generalize for any ROI
                    double referenceLat = aaps.get(0).getReferenceLat();
                    double referenceLon = aaps.get(0).getReferenceLon();

                    // Only if the reference changes, re-project
                    if (referenceLat != roiReferenceLat.get() || referenceLon != roiReferenceLon.get()) {
                        changes.getAndIncrement();
                        roiReferenceLat.set(referenceLat);
                        roiReferenceLon.set(referenceLon);
                        euclideanROI.set(Transformations.toEuclideanPlane(nonEuclideanROI, referenceLat, referenceLon));
                    }

                    List<List<double[]>> unionQueue = new ArrayList<>();

                    // For each AAP with this number of assets in sight
                    aaps.forEach(accessAreaPolygon -> {

                        List<double[]> eIntersection = new ArrayList<>();

                        try {
                            List<List<double[]>> polygonAndROI = Arrays.asList(euclideanROI.get(),
                                    Transformations.toEuclideanPlane(accessAreaPolygon.getGeoCoordinates(),
                                            referenceLat, referenceLon));
                            PolygonOperator po = new PolygonOperator();
                            Polygon intersection = po.polyIntersect(polygonAndROI, POLYGON_EPSILON);
                            eIntersection = intersection.getRegions().isEmpty() ? eIntersection : intersection.getRegions().get(0);

                        } catch (RuntimeException e) {
                            Log.error("Error trying to intersect the following polygon: ");
                            accessAreaPolygon.getGeoCoordinates().forEach(c -> Log.error(c[0] + "," + c[1]));
                            Log.error("### WITH ###");
                            nonEuclideanROI.forEach(c -> Log.error(c[0] + "," + c[1]));
                            e.printStackTrace();
                        }

                        if (eIntersection.size() >= 3) {
                            unionQueue.add(eIntersection);
                        } else if (eIntersection.size() != 0) {
                            Log.warn("Intersection with less than 3 points");
                        }

                        List<double[]> neIntersection = Transformations.toNonEuclideanPlane(eIntersection, referenceLat, referenceLon);

                        // Surface for K = 1 (intersection seen by 1 GW):
                        if (k == 1) {
                            timedMetricsRecord.addMetric(accessAreaPolygon.getGwsInSight().get(0), geographic.computeNonEuclideanSurface(neIntersection));
                        }

                        AccessAreaPolygon intersectionAccessAreaPolygon = new AccessAreaPolygon(timeElapsed, k, accessAreaPolygon.getGwsInSight(), neIntersection,
                                neIntersection.stream().map(pair -> pair[0]).collect(Collectors.toList()),
                                neIntersection.stream().map(pair -> pair[1]).collect(Collectors.toList()));
                        roiIntersections.add(intersectionAccessAreaPolygon);

                    });

                    if (!unionQueue.isEmpty()) {

                        try {
                            PolygonOperator po = new PolygonOperator();
                            Polygon union = po.polyUnion(unionQueue, POLYGON_EPSILON);
                            union = po.polyUnion(union.getRegions(), POLYGON_EPSILON);  // Second union if some polygons got clipped out
                            union.getRegions().forEach(region -> {
                                List<double[]> neIntersection = Transformations.toNonEuclideanPlane(region, referenceLat, referenceLon);
                                surfaceValues[k - 1] = surfaceValues[k - 1] + geographic.computeNonEuclideanSurface(neIntersection);
                                AccessAreaPolygon unionAccessAreaPolygon = new AccessAreaPolygon(timeElapsed, k, null, neIntersection,
                                        neIntersection.stream().map(pair -> pair[0]).collect(Collectors.toList()),
                                        neIntersection.stream().map(pair -> pair[1]).collect(Collectors.toList()));
                                roiUnions.add(unionAccessAreaPolygon);
                            });
                        } catch (NullPointerException e) {
                            Log.error("Regions empty?: " + unionQueue.isEmpty());
                            Log.error("Regions size?: " + unionQueue.size());
                            e.printStackTrace();
                        }

                    }
                });

                StringBuilder sb = new StringBuilder(timeElapsed + "");

                for (double surface : surfaceValues) {
                    sb.append(",");
                    double percentage = (surface / roiSurface) * 100D; // Math.round(((surface / roiSurface) * 100.00000) * 100000d) / 100000d;
                    sb.append(percentage);
                }

                timeSeriesData.add(timedMetricsRecord);
                statistics.add(sb.toString());

            }
        }

        Log.debug("Changes: " + changes);

        if (SAVE_GEOGRAPHIC && SAVE_SNAPSHOT) {
            saveAAPsAt(roiIntersections, "snapshot_aaps_intersection", SNAPSHOT);
            saveAAPsAt(roiUnions, "snapshot_aaps_union", SNAPSHOT);
        }

        reports.saveAsCSV(statistics, "coverage_" + (int) (satelliteList.get(0).getElements().getSemiMajorAxis() / 1000.0));
        reports.saveAsJSON(timeSeriesData, "surface_metrics");

    }

    private Map<Integer, List<AccessAreaPolygon>> mapByNOfAssets(List<AccessAreaPolygon> AccessAreaPolygons, long timeElapsed) {

        Map<Integer, List<AccessAreaPolygon>> byAssetsInSight = new LinkedHashMap<>(MAX_SUBSET_SIZE);
        AccessAreaPolygons.stream().filter(AccessAreaPolygon -> AccessAreaPolygon.getDate() == timeElapsed).forEach(AccessAreaPolygon -> {
            int nAssets = AccessAreaPolygon.getnOfGwsInSight();
            // accumulatedAreas.putIfAbsent(nAssets, surfaceInKm2);
            if (byAssetsInSight.containsKey(nAssets)) {
                byAssetsInSight.get(nAssets).add(AccessAreaPolygon);
            } else {
                List<AccessAreaPolygon> accessAreaPolygonList = new ArrayList<>();
                accessAreaPolygonList.add(AccessAreaPolygon);
                byAssetsInSight.put(nAssets, accessAreaPolygonList);
            }
        });

        return byAssetsInSight;

    }

    public void analyzeSurfaceCoverage(List<AccessAreaPolygon> AccessAreaPolygons) {

        statistics.clear();

        AbsoluteDate endDate = Utils.stamp2AD(END_DATE);
        AbsoluteDate pointerDate = Utils.stamp2AD(START_DATE);

        // Accumulated areas by number of satellites in visibility is stored in this array (idx = number of sats, value = area) // FIXME remove eventually
        Map<Integer, Double> accumulatedAreas = new HashMap<>(MAX_SUBSET_SIZE);

        while (pointerDate.compareTo(endDate) <= 0) {   // TODO while/filter can be improved

            accumulatedAreas.clear();

            long timeSinceStart = TimeUtils.stamp2unix(pointerDate.toString()) - TimeUtils.stamp2unix(START_DATE);
            AccessAreaPolygons.stream().filter(AccessAreaPolygon -> AccessAreaPolygon.getDate() == timeSinceStart).forEach(AccessAreaPolygon -> {

                double surfaceInKm2 = geographic.computeNonEuclideanSurface(AccessAreaPolygon.getGeoCoordinates()) * 1E-6; // AAP.getSurfaceInKm2(); //
                int nAssets = AccessAreaPolygon.getnOfGwsInSight();
                // accumulatedAreas.putIfAbsent(nAssets, surfaceInKm2);
                if (accumulatedAreas.containsKey(nAssets)) {
                    accumulatedAreas.put(nAssets, accumulatedAreas.get(nAssets) + surfaceInKm2);
                } else {
                    accumulatedAreas.put(nAssets, surfaceInKm2);
                }
                AccessAreaPolygon.setSurfaceInKm2(surfaceInKm2);

            });

            // Save the results // FIXME Define a result object
            statistics.add(stringifyResults(pointerDate, accumulatedAreas));

            // Advance to the next time step
            pointerDate = pointerDate.shiftedBy(TIME_STEP);

        }

        reports.saveAsCSV(statistics, "stats");

    }

    private void saveAAPs(List<AccessAreaPolygon> AccessAreaPolygons, String fileName) {

        for (int nOfGw = 1; nOfGw <= MAX_SUBSET_SIZE; nOfGw++) {
            int finalNOfGw = nOfGw;

            try {
                reports.saveAsJSON(AccessAreaPolygons.stream()
                        .filter(AccessAreaPolygon -> AccessAreaPolygon.getnOfGwsInSight() == finalNOfGw)
                        .collect(Collectors.toList()), fileName + "_" + nOfGw);
            } catch (IllegalArgumentException e) {
                AccessAreaPolygons.stream()
                        .filter(AccessAreaPolygon -> AccessAreaPolygon.getnOfGwsInSight() == finalNOfGw)
                        .collect(Collectors.toList()).forEach(accessAreaPolygon ->
                                Log.error(accessAreaPolygon.getReferenceLat() + "," + accessAreaPolygon.getReferenceLon() + "," + accessAreaPolygon.getSurfaceInKm2()));
            }
        }

    }

    private void saveAAPsAt(List<AccessAreaPolygon> AccessAreaPolygons, String fileName, long time) {
        reports.saveAsJSON(AccessAreaPolygons.stream()
                .filter(AccessAreaPolygon -> AccessAreaPolygon.getDate() == time)
                .collect(Collectors.toList()), fileName);
    }

    // TODO REPLACE WITH POST-ANALYSIS
    private String stringifyResults(AbsoluteDate pointerDate, Map<Integer, Double> accumulatedAreas) {

        StringBuilder sb = new StringBuilder();
        sb.append(TimeUtils.stamp2unix(pointerDate.toString()) - UNIX_START_DATE);

        for (Integer key : accumulatedAreas.keySet()) {
            sb.append(",");
            sb.append(accumulatedAreas.get(key));
        }

        return sb.toString();

    }

    /**
     * This method takes the satellite list, propagates orbits to the specified date and computes the corresponding
     * access area or FOV polygon for each one.
     *
     * @param satelliteList a list of Satellite objects
     * @param date          an AbsoluteDate object
     * @return a List of Regions
     **/
    private List<AccessRegion> computeAccessRegionsAt(List<Satellite> satelliteList, Simulation simulation, AbsoluteDate date, long timeElapsed) {

        List<AccessRegion> AccessRegionList = new ArrayList<>();

        for (Satellite satellite : satelliteList) {

            Ephemeris eph;

            if (PROPAGATE_INTERNALLY) {
                simulation.setSatellite(satellite);
                eph = simulation.computeFixedEphemerisKm(date);
                eph.setPos(eph.getPosX(), eph.getPosY(), eph.getPosZ());
                eph.setSSP(Math.toDegrees(eph.getLatitude()), Math.toDegrees(eph.getLongitude()), eph.getHeight());
            } else {
                eph = constellation.get(satellite.getId()).get(timeElapsed);
            }

            double x = 0, y = 0, z = 0;

            try {
                x = eph.getPosX();
                y = eph.getPosY();
                z = eph.getPosZ();
            } catch (NullPointerException e) {
                Log.error("time: " + timeElapsed);
                Log.error("sat id: " + satellite.getId());
                e.printStackTrace();
            }

            double lambdaMax = geographic.getLambdaMax(x, y, z, VISIBILITY_THRESHOLD);

            List<double[]> poly = OblateAccessRegion.drawLLAConic(x, y, z, VISIBILITY_THRESHOLD, 1E-4, POLYGON_SEGMENTS);
//            List<double[]> poly = OblateFOV.drawLLAConic(x, y, z, 69, POLYGON_SEGMENTS);
//             List<double[]> poly = geo.drawSphericalAAP(lambdaMax, eph.getLatitude(), eph.getLongitude(), POLYGON_SEGMENTS);

            AccessRegion accessRegion = new AccessRegion();
            accessRegion.setLambdaMax(lambdaMax);
            accessRegion.setSatLLA(eph.getLatitude(), eph.getLongitude(), eph.getHeight());
            accessRegion.setSatPos(x, y, z);
            accessRegion.setPolygonCoordinates(poly);
            AccessRegionList.add(accessRegion);

        }

        return AccessRegionList;

    }

    /**
     * This method checks the distance between each and all pairs of coordinates corresponding to the satellites SSPs,
     * and returns whether the intersection is empty or not
     *
     * @param assetsToCheck    A List of Integer depicting the indexes of the satellites that are being checked
     * @param AccessRegionList A List of Region objects that will provide the SSP coordinates for the assets being checked
     * @return boolean whether the intersection is empty or not
     **/
    private boolean checkDistances(List<Integer> assetsToCheck, List<AccessRegion> AccessRegionList) {

        // First we need to generate the combination list for the pair of assets that need checking
        Combination combination = new Combination();

        // TODO: Implement a quick discard feature that excludes regions totally unable to intersect

        List<List<Integer>> pairsToCheck = combination.computeCombinations(assetsToCheck, 2);

        for (List<Integer> pairToCheck : pairsToCheck) {

            int r1Idx = pairToCheck.get(0);
            int r2Idx = pairToCheck.get(1);

            double lambda1 = AccessRegionList.get(r1Idx).getLambdaMax();
            double lambda2 = AccessRegionList.get(r2Idx).getLambdaMax();
            double distance = geographic.computeGeodesic(AccessRegionList.get(r1Idx), AccessRegionList.get(r2Idx));

            if (distance >= (lambda1 + lambda2) * (1 + LAMBDA_EXCLUSION / 100.0)) {
                return false;
            }
        }
        return true;
    }

}