package constellation.tools;

import com.menecats.polybool.models.Polygon;
import constellation.tools.geometry.AccessAreaPolygon;
import constellation.tools.geometry.AccessRegion;
import constellation.tools.geometry.GeographicTools;
import constellation.tools.geometry.OblateAccessRegion;
import constellation.tools.math.Combination;
import constellation.tools.math.TimedMetricsRecord;
import constellation.tools.math.Transformations;
import constellation.tools.operations.PolygonOperator;
import constellation.tools.utilities.AppConfig;
import constellation.tools.utilities.ConstellationHash;
import constellation.tools.utilities.FileUtils;
import constellation.tools.utilities.TimeUtils;
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

import java.io.File;
import java.util.*;
import java.util.concurrent.CopyOnWriteArrayList;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.AtomicReference;

/**
 * Dynamic Constellation Coverage Computer (D3CO)
 **/
public class ConstellationCoverageComputer {

    private List<Satellite> satelliteList;
    private final List<String> statistics = new ArrayList<>();
    private List<Map<Long, Ephemeris>> constellation;
    private List<double[]> nonEuclideanROI;

    private List<TimedMetricsRecord> constellationCoverageTimeSeries;

    private static AppConfig appConfig;
    private Simulation simulation;
    private FileUtils fileUtils;
    private GeographicTools geographicTools;
    private PolygonOperator polygonOperator;
    private String simulationHash;

    /* Metrics */
    int avoidedDueToNotIntersectingCount = 0;
    int polyOperationsPerformedCount = 0;
    int polyOperationsThatResultedInEmptyIntersection = 0;

    // TODO: Add ban list

    /**
     * Default constructor
     **/
    public ConstellationCoverageComputer() {
        loadConfigurations();
        satelliteList = Utils.satellitesFromFile(appConfig.satellitesFile());
    }

    /**
     * Accepts configurations (either config file or JSON)
     **/
    public ConstellationCoverageComputer(String configurations) {
        loadConfigurations(configurations);
        satelliteList = Utils.satellitesFromFile(appConfig.satellitesFile());
    }

    /**
     * Attempts to locate and load configuration.properties on root
     **/
    public void loadConfigurations() {
        loadConfigurations("config.properties");
    }

    public void loadConfigurations(String configurationsFilePath) {
        try {
            appConfig = new AppConfig(configurationsFilePath);
            fileUtils = new FileUtils(appConfig.outputPath());
            simulation = new Simulation(appConfig.orekitPath());
            simulation.setParams(appConfig.startDate(), appConfig.endDate(), appConfig.timeStep(), appConfig.visibilityThreshold());
            simulation.setInertialFrame(FramesFactory.getEME2000());
            simulation.setFixedFrame(FramesFactory.getITRF(IERSConventions.IERS_2010, true));
            polygonOperator = new PolygonOperator(appConfig.polygonEpsilon());

            if (appConfig.useRoiFile()) {
                nonEuclideanROI = FileUtils.file2DoubleList(appConfig.roiPath());
            }
        } catch (ConfigurationException e) {
            Log.error("Error trying to load configurations. Exiting. " + e.getMessage());
            System.exit(100);
        } catch (Exception e) {
            Log.error("Error trying to instance the simulation. Exiting. ");
            e.printStackTrace();
            System.exit(101);
        }

        geographicTools = new GeographicTools();
    }

    public void run() {

        setSimulationHash(ConstellationHash.hash2(satelliteList));
        Log.info("Simulation hash: " + this.simulationHash);

        if (!appConfig.propagateInternally()) {
            Log.info("Reading ephemeris from directory: " + appConfig.positionsPath());
            constellation = fileUtils.positionsFromPath(appConfig.positionsPath(), satelliteList.size());
        }

        List<AccessAreaPolygon> nonEuclideanAccessAreaPolygons = propagateAccessAreaPolygons();
        analyzeROICoverage(nonEuclideanAccessAreaPolygons);

    }

    public void computeEarthCoverage() {

        setSimulationHash(ConstellationHash.hash2(satelliteList));
        Log.info("Simulation hash: " + this.simulationHash);

        if (!appConfig.propagateInternally()) {
            Log.debug("Reading ephemeris from directory: " + appConfig.positionsPath());
            constellation = fileUtils.positionsFromPath(appConfig.positionsPath(), satelliteList.size());
        }

        List<AccessAreaPolygon> nonEuclideanAccessAreaPolygons = propagateAccessAreaPolygons();
        analyzeConstellationCoverage(nonEuclideanAccessAreaPolygons);

    }

    private List<AccessAreaPolygon> propagateAccessAreaPolygons() {

        AbsoluteDate startDate = Utils.stamp2AD(appConfig.startDate());
        AbsoluteDate endDate = Utils.stamp2AD(appConfig.endDate(), DataContext.getDefault().getTimeScales().getUTC());

        // We compute a Utility "List of Lists", containing all possible overlapping combinations between regions.
        Combination comb = new Combination(satelliteList.size(), Math.min(appConfig.maxSubsetSize(), satelliteList.size()));
        final List<List<Integer>> combinationsList = comb.computeCombinations();

        List<AccessAreaPolygon> nonEuclideanAccessAreaPolygons = new ArrayList<>();
        List<AccessAreaPolygon> euclideanAccessAreaPolygons = new ArrayList<>();

        for (AbsoluteDate date = startDate; date.compareTo(endDate) <= 0; date = date.shiftedBy(appConfig.timeStep())) {

            long timeElapsed = TimeUtils.stamp2unix(date.toString()) - TimeUtils.stamp2unix(appConfig.startDate());

            // Obtain the starting non-euclidean FOVs
            List<AccessRegion> nonEuclideanAccessRegions = computeAccessRegionsAt(satelliteList, simulation, date, timeElapsed);
            computeIntersections(combinationsList, nonEuclideanAccessAreaPolygons, euclideanAccessAreaPolygons, timeElapsed, nonEuclideanAccessRegions);
            Log.debug("Performed: " + polyOperationsPerformedCount + " - Avoided: " + avoidedDueToNotIntersectingCount + " - Escaped: " + polyOperationsThatResultedInEmptyIntersection);

        }

        writePolygonsToFile(nonEuclideanAccessAreaPolygons, euclideanAccessAreaPolygons);

        return nonEuclideanAccessAreaPolygons;
    }

    private void writePolygonsToFile(List<AccessAreaPolygon> nonEuclideanAccessAreaPolygons, List<AccessAreaPolygon> euclideanAccessAreaPolygons) {

        if (appConfig.saveGeographic())
            saveAAPs(nonEuclideanAccessAreaPolygons, "ne_polygons_" + simulationHash);

        if (appConfig.saveEuclidean())
            saveAAPs(euclideanAccessAreaPolygons, "e_polygons_" + simulationHash);

        if (appConfig.saveSnapshot())
            saveAAPsAt(nonEuclideanAccessAreaPolygons, "snapshot_ne_polygons_" + simulationHash, appConfig.snapshot());

        if (appConfig.saveSnapshot() && appConfig.saveEuclidean())
            saveAAPsAt(euclideanAccessAreaPolygons, "snapshot_e_polygons_" + simulationHash, appConfig.snapshot());

    }

    private void computeIntersections(List<List<Integer>> combinationsList, List<AccessAreaPolygon> nonEuclideanAccessAreaPolygons, List<AccessAreaPolygon> euclideanAccessAreaPolygons, long timeElapsed, List<AccessRegion> nonEuclideanAccessRegions) {

        for (List<Integer> combination : combinationsList) {

            List<double[]> nonEuclideanCoordinates = new ArrayList<>();
            List<double[]> euclideanCoordinates = new ArrayList<>();

            // assemble a list of the FOVs to be intersected at this time step:
            List<AccessRegion> accessRegionsPendingIntersection = new ArrayList<>();
            combination.forEach(regionIndex -> accessRegionsPendingIntersection.add(nonEuclideanAccessRegions.get(regionIndex)));

            // TODO use satellite lambda
            double referenceLat = 0;
            double referenceLon = 0;

            // If this is the immediate FOV for a single satellite
            if (combination.size() <= 1) {
                int fovIdx = combination.get(0);
                AccessRegion neAccessRegion = nonEuclideanAccessRegions.get(fovIdx);
                nonEuclideanCoordinates = nonEuclideanAccessRegions.get(fovIdx).getPolygonCoordinates();
                euclideanCoordinates = Transformations.toEuclideanPlane(neAccessRegion.getPolygonCoordinates(),
                        referenceLat, referenceLon);

            } else if (checkIntersection(combination, nonEuclideanAccessRegions) && !accessRegionsPendingIntersection.isEmpty()) {

                polyOperationsPerformedCount++;
                List<List<double[]>> intersectionQueue = new LinkedList<>();
                accessRegionsPendingIntersection.forEach(accessRegion -> intersectionQueue.add(Transformations.toEuclideanPlane(accessRegion.getPolygonCoordinates(), referenceLat, referenceLon)));

                // Obtain access polygon
                Polygon intersection = polygonOperator.polyIntersect(intersectionQueue);
                if (!intersection.getRegions().isEmpty() && intersection.getRegions().get(0).size() > 2) {
                    nonEuclideanCoordinates = Transformations.toNonEuclideanPlane(intersection.getRegions().get(0),
                            referenceLat, referenceLon);
                    if (appConfig.saveEuclidean()) {
                        euclideanCoordinates = intersection.getRegions().get(0);
                    }
                } else {
                    polyOperationsThatResultedInEmptyIntersection++;
                }

            } else {
                avoidedDueToNotIntersectingCount++;
            }

            recordPolygon(nonEuclideanAccessAreaPolygons, timeElapsed, combination, nonEuclideanCoordinates, referenceLat, referenceLon);

            if (appConfig.saveEuclidean()) {
                recordPolygon(euclideanAccessAreaPolygons, timeElapsed, combination, euclideanCoordinates, referenceLat, referenceLon);
            }

        }
    }

    private void recordPolygon(List<AccessAreaPolygon> accessAreaPolygonList,
                               long timeElapsed, List<Integer> combination, List<double[]> coordinates,
                               double referenceLat, double referenceLon) {
        // FIX: Return early for degenerate/empty polygons instead of storing them.
        // Storing empty polygons causes polyIntersect(roi, emptyPolygon) to return the
        // full ROI in some PolygonBool library versions, falsely reporting 100% coverage.
        if (coordinates.size() <= 2) {
            Log.debug("Less than 2 coordinates in polygon intersection — skipping record.");
            return;
        }
        accessAreaPolygonList.add(new AccessAreaPolygon(timeElapsed, combination.size(), combination, coordinates,
                coordinates.stream().map(pair -> pair[0]).toList(),
                coordinates.stream().map(pair -> pair[1]).toList(),
                referenceLat, referenceLon));
    }

    public void analyzeROICoverage(List<AccessAreaPolygon> accessAreaPolygons) {

        statistics.clear();

        double roiSurface = geographicTools.computeNonEuclideanSurface(nonEuclideanROI);
        Log.info("ROI Surface: " + roiSurface);

        AbsoluteDate startDate = Utils.stamp2AD(appConfig.startDate());
        AbsoluteDate endDate = Utils.stamp2AD(appConfig.endDate());
        long startTimestamp = TimeUtils.stamp2unix(appConfig.startDate());

        // FIX: Use CopyOnWriteArrayList — these lists are written from the inner stream
        // and read later for snapshot saving. Plain ArrayList is not safe here.
        List<AccessAreaPolygon> roiIntersections = new CopyOnWriteArrayList<>();
        List<AccessAreaPolygon> roiUnions = new CopyOnWriteArrayList<>();
        this.constellationCoverageTimeSeries = new ArrayList<>();
        AtomicInteger changes = new AtomicInteger();

        Log.info("Analyzing ROI coverage");
        for (AbsoluteDate t = startDate; t.compareTo(endDate) <= 0; t = t.shiftedBy(appConfig.timeStep())) {

            long timeElapsed = TimeUtils.stamp2unix(t.toString()) - startTimestamp;

            Map<Integer, CopyOnWriteArrayList<AccessAreaPolygon>> byAssetsInSight = mapByNOfAssets(accessAreaPolygons, timeElapsed);

            // FIX: Guard against missing k=1 entry. If no single-satellite AAP exists at this
            // timestep (all satellites below horizon), skip rather than NPE on .get(1).get(0).
            if (!byAssetsInSight.containsKey(1) || byAssetsInSight.get(1).isEmpty()) {
                Log.warn("No single-satellite AAPs at t=" + timeElapsed + "; skipping timestep.");
                constellationCoverageTimeSeries.add(new TimedMetricsRecord(timeElapsed, satelliteList.size()));
                statistics.add(timeElapsed + ",".repeat(satelliteList.size()));
                continue;
            }

            // FIX: surfaceValues is written from the byAssetsInSight.forEach lambda (sequential),
            // so a plain double[] is safe as long as the inner stream is also sequential (see below).
            double[] surfaceValues = new double[satelliteList.size()];

            AtomicReference<Double> roiReferenceLat = new AtomicReference<>(byAssetsInSight.get(1).get(0).getReferenceLat());
            AtomicReference<Double> roiReferenceLon = new AtomicReference<>(byAssetsInSight.get(1).get(0).getReferenceLon());
            AtomicReference<List<double[]>> euclideanROI = new AtomicReference<>(
                    Transformations.toEuclideanPlane(nonEuclideanROI, roiReferenceLat.get(), roiReferenceLon.get()));

            TimedMetricsRecord timedMetricsRecord = new TimedMetricsRecord(timeElapsed, satelliteList.size());

            byAssetsInSight.forEach((k, aaps) -> {

                // FIX: Validate k is within array bounds before using as index.
                // k comes from getnOfGwsInSight() which reflects combination.size(); a config
                // mismatch (maxSubsetSize > satelliteList.size()) can produce k > array length.
                if (k < 1 || k > satelliteList.size()) {
                    Log.warn("Unexpected k=" + k + " (satelliteList.size()=" + satelliteList.size() + "); skipping group.");
                    return;
                }

                double referenceLat = aaps.get(0).getReferenceLat();
                double referenceLon = aaps.get(0).getReferenceLon();

                if (referenceLat != roiReferenceLat.get() || referenceLon != roiReferenceLon.get()) {
                    changes.getAndIncrement();
                    roiReferenceLat.set(referenceLat);
                    roiReferenceLon.set(referenceLon);
                    euclideanROI.set(Transformations.toEuclideanPlane(nonEuclideanROI, referenceLat, referenceLon));
                }

                // FIX: Changed to plain ArrayList — safe because the stream below is sequential.
                // The original parallelStream() caused concurrent ArrayList.add() calls, which
                // corrupts internal state and throws ArrayIndexOutOfBoundsException.
                List<List<double[]>> unionQueue = new ArrayList<>();

                // FIX: Changed parallelStream() → stream() (sequential).
                // parallelStream() caused concurrent mutation of: unionQueue (ArrayList),
                // surfaceValues (double[]), timedMetricsRecord, and roiIntersections (ArrayList).
                // None of these are safe for concurrent writes.
                aaps.stream().forEach(accessAreaPolygon -> {

                    try {
                        // FIX: Skip AAPs with degenerate/empty geometry before calling polyIntersect.
                        // An empty polygon fed into polyIntersect may be treated as the universe by
                        // the PolygonBool library, returning the full ROI and falsely reporting 100%.
                        List<double[]> aapGeoCoords = accessAreaPolygon.getGeoCoordinates();
                        if (aapGeoCoords == null || aapGeoCoords.size() < 3) {
                            Log.debug("Skipping AAP with < 3 geo coordinates at k=" + k + ", t=" + timeElapsed);
                            return;
                        }

                        // FIX: Validate the Euclidean projection also has enough points before
                        // passing to polyIntersect. toEuclideanPlane mirrors input length, but
                        // guard explicitly in case upstream data is malformed.
                        List<double[]> euclideanAAP = Transformations.toEuclideanPlane(aapGeoCoords, referenceLat, referenceLon);
                        if (euclideanAAP.size() < 3) {
                            Log.debug("Projected AAP has < 3 points at k=" + k + ", t=" + timeElapsed + "; skipping.");
                            return;
                        }

                        List<List<double[]>> polygonAndROI = Arrays.asList(
                                euclideanROI.get(),
                                euclideanAAP
                        );
                        Polygon intersection = polygonOperator.polyIntersect(polygonAndROI);

                        // FIX: Multipart intersection support — previously only .get(0) was used,
                        // silently dropping any additional disjoint regions. Now all parts are processed.
                        List<List<double[]>> allRegions = intersection.getRegions();

                        if (allRegions.isEmpty()) {
                            // No intersection with ROI — nothing to record for this AAP.
                            return;
                        }

                        if (allRegions.size() > 1) {
                            Log.debug("Multipart intersection: " + allRegions.size() + " parts at k=" + k + ", t=" + timeElapsed);
                        }

                        for (List<double[]> eIntersectionPart : allRegions) {

                            // FIX: Skip degenerate regions — fewer than 3 points cannot form a polygon.
                            if (eIntersectionPart.size() < 3) {
                                Log.warn("Intersection part has < 3 points; skipping.");
                                continue;
                            }

                            List<double[]> neIntersectionPart = Transformations.toNonEuclideanPlane(eIntersectionPart, referenceLat, referenceLon);

                            // Per-satellite surface contribution (k=1 only)
                            if (k == 1) {
                                // FIX: Guard gwsInSight before indexing to avoid IndexOutOfBoundsException.
                                List<Integer> gwsInSight = accessAreaPolygon.getGwsInSight();
                                if (gwsInSight != null && !gwsInSight.isEmpty()) {
                                    timedMetricsRecord.addMetric(gwsInSight.get(0),
                                            geographicTools.computeNonEuclideanSurface(neIntersectionPart));
                                }
                            }

                            // Each disjoint part is queued individually for the union pass
                            unionQueue.add(eIntersectionPart);

                            AccessAreaPolygon intersectionAAP = new AccessAreaPolygon(
                                    timeElapsed, k, accessAreaPolygon.getGwsInSight(), neIntersectionPart,
                                    neIntersectionPart.stream().map(pair -> pair[0]).toList(),
                                    neIntersectionPart.stream().map(pair -> pair[1]).toList());
                            roiIntersections.add(intersectionAAP);
                        }

                    } catch (RuntimeException e) {
                        Log.error("Error intersecting polygon at k=" + k + ", t=" + timeElapsed + ": " + e.getMessage());
                        accessAreaPolygon.getGeoCoordinates().forEach(c -> Log.error(c[0] + "," + c[1]));
                        Log.error("### WITH ROI ###");
                        nonEuclideanROI.forEach(c -> Log.error(c[0] + "," + c[1]));
                    }
                });

                if (!unionQueue.isEmpty()) {
                    try {
                        Polygon union = polygonOperator.polyUnion(unionQueue);
                        union = polygonOperator.polyUnion(union.getRegions()); // Second pass for clipped polygons

                        for (List<double[]> region : union.getRegions()) {
                            // FIX: Skip degenerate union regions.
                            if (region.size() < 3) {
                                Log.warn("Union region has < 3 points at k=" + k + "; skipping.");
                                continue;
                            }
                            List<double[]> neUnionRegion = Transformations.toNonEuclideanPlane(region, referenceLat, referenceLon);
                            // k already validated above — bounds-safe.
                            surfaceValues[k - 1] += geographicTools.computeNonEuclideanSurface(neUnionRegion);

                            AccessAreaPolygon unionAAP = new AccessAreaPolygon(
                                    timeElapsed, k, null, neUnionRegion,
                                    neUnionRegion.stream().map(pair -> pair[0]).toList(),
                                    neUnionRegion.stream().map(pair -> pair[1]).toList());
                            roiUnions.add(unionAAP);
                        }
                    } catch (NullPointerException e) {
                        Log.error("Union failed at k=" + k + ", t=" + timeElapsed
                                + "; unionQueue.size()=" + unionQueue.size() + ": " + e.getMessage());
                    }
                }
            });

            StringBuilder sb = new StringBuilder(String.valueOf(timeElapsed));
            for (double surface : surfaceValues) {
                sb.append(",");
                sb.append((surface / roiSurface) * 100D);
            }

            timedMetricsRecord.scale(1.0 / roiSurface);
            constellationCoverageTimeSeries.add(timedMetricsRecord);
            statistics.add(sb.toString());
        }

        Log.info("Ending ROI coverage analysis");
        Log.debug("Reprojection changes: " + changes);

        if (appConfig.saveGeographic() && appConfig.saveSnapshot()) {
            saveAAPsAt(roiIntersections, "snapshot_aaps_intersection_" + simulationHash, appConfig.snapshot());
            saveAAPsAt(roiUnions, "snapshot_aaps_union_" + simulationHash, appConfig.snapshot());
        }

        fileUtils.saveAsCSV(statistics, "coverage_" + simulationHash);
        fileUtils.saveAsJSON(constellationCoverageTimeSeries, "surface_metrics_" + simulationHash);
    }

    public void analyzeConstellationCoverage(List<AccessAreaPolygon> accessAreaPolygons) {

        statistics.clear();

        // Timekeeping
        AbsoluteDate startDate = Utils.stamp2AD(appConfig.startDate());
        AbsoluteDate endDate = Utils.stamp2AD(appConfig.endDate());
        long startTimestamp = TimeUtils.stamp2unix(appConfig.startDate());

        List<TimedMetricsRecord> timeSeriesData = new ArrayList<>();

        for (AbsoluteDate t = startDate; t.compareTo(endDate) <= 0; t = t.shiftedBy(appConfig.timeStep())) {

            long timeElapsed = TimeUtils.stamp2unix(t.toString()) - startTimestamp;

            // Group regions by number of satellites on sight, for this particular time step
            Map<Integer, CopyOnWriteArrayList<AccessAreaPolygon>> byAssetsInSight = mapByNOfAssets(accessAreaPolygons, timeElapsed);

            double[] surfaceValues = new double[satelliteList.size()];
            Arrays.fill(surfaceValues, 0.0);
            TimedMetricsRecord timedMetricsRecord = new TimedMetricsRecord(timeElapsed, satelliteList.size());

            byAssetsInSight.forEach((k, aaps) -> aaps.forEach(accessAreaPolygon -> {
                surfaceValues[k - 1] += geographicTools.computeNonEuclideanSurface(accessAreaPolygon.getGeoCoordinates());
                timedMetricsRecord.accumulateMetric(k - 1, geographicTools.computeNonEuclideanSurface(accessAreaPolygon.getGeoCoordinates()));
            }));

            StringBuilder sb = new StringBuilder(timeElapsed + "");
            for (double surface : surfaceValues) {
                sb.append(",");
                sb.append(surface);
            }

            timedMetricsRecord.scale(10.0E-6); // To m2
            timeSeriesData.add(timedMetricsRecord);
            statistics.add(sb.toString());

        }
        Log.info("Ending Constellation Coverage Analysis");

//        if (appConfig.saveGeographic() && appConfig.saveSnapshot()) {
//            saveAAPsAt(roiIntersections, "snapshot_aaps_intersection_" + simulationHash, appConfig.snapshot());
//            saveAAPsAt(roiUnions, "snapshot_aaps_union_" + simulationHash, appConfig.snapshot());
//        }

        // fileUtils.saveAsCSV(statistics, "constellation_coverage_" + simulationHash);
        fileUtils.saveAsCSV(timeSeriesData.stream().map(TimedMetricsRecord::toString).toList(), "constellation_coverage_" + simulationHash);
        fileUtils.saveAsJSON(timeSeriesData, "constellation_surface_metrics_" + simulationHash);

    }

    private Map<Integer, CopyOnWriteArrayList<AccessAreaPolygon>> mapByNOfAssets(List<AccessAreaPolygon> accessAreaPolygons, long timeElapsed) {

        Map<Integer, CopyOnWriteArrayList<AccessAreaPolygon>> byAssetsInSight = new LinkedHashMap<>(Math.min(appConfig.maxSubsetSize(), satelliteList.size()));
        accessAreaPolygons.stream().filter(aap -> aap.getDate() == timeElapsed).forEach(aap -> {
            int nAssets = aap.getnOfGwsInSight();
            if (byAssetsInSight.containsKey(nAssets)) {
                byAssetsInSight.get(nAssets).add(aap);
            } else {
                CopyOnWriteArrayList<AccessAreaPolygon> accessAreaPolygonList = new CopyOnWriteArrayList<>();
                accessAreaPolygonList.add(aap);
                byAssetsInSight.put(nAssets, accessAreaPolygonList);
            }
        });

        return byAssetsInSight;

    }

    private void saveAAPs(List<AccessAreaPolygon> accessAreaPolygons, String fileName) {

        for (int nOfGw = 1; nOfGw <= Math.min(appConfig.maxSubsetSize(), satelliteList.size()); nOfGw++) {
            int finalNOfGw = nOfGw;

            try {
                fileUtils.saveAsJSON(accessAreaPolygons.stream()
                        .filter(aap -> aap.getnOfGwsInSight() == finalNOfGw).toList(), fileName + "_" + nOfGw);
            } catch (IllegalArgumentException e) {
                accessAreaPolygons.stream()
                        .filter(aap -> aap.getnOfGwsInSight() == finalNOfGw)
                        .toList().forEach(accessAreaPolygon ->
                                Log.error(accessAreaPolygon.getReferenceLat()
                                        + "," + accessAreaPolygon.getReferenceLon()
                                        + "," + accessAreaPolygon.getSurfaceInKm2()));
            }
        }

    }

    private void saveAAPsAt(List<AccessAreaPolygon> accessAreaPolygons, String fileName, long time) {
        fileUtils.saveAsJSON(accessAreaPolygons.stream()
                .filter(aap -> aap.getDate() == time).toList(), fileName);
    }

    /**
     * This method takes the satellite list, propagates orbits to the specified date and computes the corresponding
     * Oblate Access Region for each one.
     *
     * @param satelliteList a list of Satellite objects
     * @param date          an AbsoluteDate object
     * @return a List of Regions
     **/
    private List<AccessRegion> computeAccessRegionsAt(List<Satellite> satelliteList, Simulation simulation, AbsoluteDate date, long timeElapsed) {

        List<AccessRegion> accessRegions = new ArrayList<>();

        for (Satellite satellite : satelliteList) {

            Ephemeris eph;

            if (appConfig.propagateInternally()) {
                simulation.setSatellite(satellite);
                eph = simulation.computeFixedEphemerisKm(date);
                eph.setPos(eph.getPosX(), eph.getPosY(), eph.getPosZ());
                eph.setSSP(Math.toDegrees(eph.getLatitude()), Math.toDegrees(eph.getLongitude()), eph.getHeight());
            } else {
                eph = constellation.get(satellite.getId()).get(timeElapsed);
            }

            double x = 0;
            double y = 0;
            double z = 0;

            try {
                x = eph.getPosX();
                y = eph.getPosY();
                z = eph.getPosZ();
            } catch (NullPointerException e) {
                Log.error("time: " + timeElapsed);
                Log.error("sat id: " + satellite.getId());
                Log.error(e.getLocalizedMessage());
            }

            double lambdaMax = geographicTools.getLambdaMax(x, y, z, appConfig.visibilityThreshold());

            List<double[]> poly = OblateAccessRegion.drawLLAConic(x, y, z, appConfig.visibilityThreshold(), 1E-4, appConfig.polygonSegments());

            AccessRegion accessRegion = new AccessRegion();
            accessRegion.setLambdaMax(lambdaMax);
            accessRegion.setSatLLA(eph.getLatitude(), eph.getLongitude(), eph.getHeight());
            accessRegion.setSatPos(x, y, z);
            accessRegion.setPolygonCoordinates(poly);
            accessRegions.add(accessRegion);

        }

        return accessRegions;

    }

    /**
     * This method checks the distance between each and all pairs of coordinates corresponding to the satellites SSPs,
     * and returns whether the intersection is empty or not based on the SSP and Earth's Max Central Angle Lambda
     *
     * @param assetsToCheck A List of Integer depicting the indexes of the satellites that are being checked
     * @param accessRegions A List of Region objects that will provide the SSP coordinates for the assets being checked
     * @return boolean whether the intersection is empty or not
     **/
    private boolean checkIntersection(List<Integer> assetsToCheck, List<AccessRegion> accessRegions) {

        // First we need to generate the combination list for the pair of assets that need checking
        Combination combination = new Combination();

        // TODO: Implement a quick discard feature that excludes regions totally unable to intersect

        List<List<Integer>> pairsToCheck = combination.computeCombinations(assetsToCheck, 2);

        for (List<Integer> pairToCheck : pairsToCheck) {

            int r1Idx = pairToCheck.get(0);
            int r2Idx = pairToCheck.get(1);

            double lambda1 = accessRegions.get(r1Idx).getLambdaMax();
            double lambda2 = accessRegions.get(r2Idx).getLambdaMax();
            double distance = geographicTools.computeGeodesic(accessRegions.get(r1Idx), accessRegions.get(r2Idx));

            if (distance >= (lambda1 + lambda2) * (1 + appConfig.lambdaExclusion() / 100.0)) {
                return false;
            }
        }
        return true;
    }

    public void setOutputPath(String outputPath) {

        File directory = new File(outputPath);

        if (!directory.exists()) {
            if (directory.mkdirs()) {
                Log.info("Output directory didn't exist, so it was created.");
            } else {
                Log.error("Failed to create the output directory.");
            }
        }

        fileUtils = new FileUtils(outputPath);

    }

    public void setROI(List<double[]> roiPolygonInGeoCoordinates) {
        this.nonEuclideanROI = roiPolygonInGeoCoordinates;
    }

    public List<Satellite> getSatelliteList() {
        return satelliteList;
    }

    public void setSatelliteList(List<Satellite> satelliteList) {
        this.satelliteList = satelliteList;
    }

    private void setSimulationHash(String simulationHash) {
        this.simulationHash = simulationHash;
    }

    public String getSimulationHash() {
        return this.simulationHash;
    }

    public List<TimedMetricsRecord> getConstellationCoverageTimeSeries() {
        return this.constellationCoverageTimeSeries;
    }

}