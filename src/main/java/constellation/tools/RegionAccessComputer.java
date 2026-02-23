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
 * Region Access Computer uses the same set of tools as the ConstellationCoverageComputer, but deals with the access metrics
 * for the satellite itself, as an object, rather than its FOV.
 **/
public class RegionAccessComputer {

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
    public RegionAccessComputer() {
        loadConfigurations();
        satelliteList = Utils.satellitesFromFile(appConfig.satellitesFile());
    }

    /**
     * Accepts configurations (either config file or JSON)
     **/
    public RegionAccessComputer(String configurations) {
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
            Log.error(e.getLocalizedMessage());
            System.exit(101);
        }

        geographicTools = new GeographicTools();
    }

    public void computeROIAccess() {

        setSimulationHash(ConstellationHash.hash2(satelliteList));
        Log.info("Simulation hash: " + this.simulationHash);

        if (!appConfig.propagateInternally()) {
            Log.debug("Reading ephemeris from directory: " + appConfig.positionsPath());
            constellation = fileUtils.positionsFromPath(appConfig.positionsPath(), satelliteList.size());
        }

        List<TimedMetricsRecord> roiAccessStatistics = propagateOrbitsAndObtainROIIntrusion();
        fileUtils.saveAsJSON(roiAccessStatistics, "roi_access_metrics_" + simulationHash);

    }

    public void computeROMetrics() {

        setSimulationHash(ConstellationHash.hash2(satelliteList));
        Log.info("Simulation hash: " + this.simulationHash);

        if (!appConfig.propagateInternally()) {
            Log.debug("Reading ephemeris from directory: " + appConfig.positionsPath());
            constellation = fileUtils.positionsFromPath(appConfig.positionsPath(), satelliteList.size());
        }

        List<TimedMetricsRecord> roiAccessStatistics = propagateOrbitsAndObtainROIMetrics();
        fileUtils.saveAsJSON(roiAccessStatistics, "roi_access_metrics_" + simulationHash);

    }

    // TODO: Now, I'm computing everything for time's sake and not wanting to mess up anything. Access Area Polygons include
    // the SSP due to Santi in the past actions.
    private List<TimedMetricsRecord> propagateOrbitsAndObtainROIIntrusion() {

        AbsoluteDate endDate = Utils.stamp2AD(appConfig.endDate(), DataContext.getDefault().getTimeScales().getUTC());

        List<TimedMetricsRecord> roiAccessStatistics = new ArrayList<>();

        for (AbsoluteDate date = Utils.stamp2AD(appConfig.startDate()); date.compareTo(endDate) <= 0; date = date.shiftedBy(appConfig.timeStep())) {

            long timeElapsed = TimeUtils.stamp2unix(date.toString()) - TimeUtils.stamp2unix(appConfig.startDate());

            // Obtain the starting non-euclidean FOVs
            List<AccessRegion> nonEuclideanAccessRegions = computeAccessRegionsAt(satelliteList, simulation, date, timeElapsed);

            // Check satellites within the ROI
            computeSatellitesWithinROI(roiAccessStatistics, nonEuclideanAccessRegions, timeElapsed);

        }

        return roiAccessStatistics;
    }

    private List<TimedMetricsRecord> propagateOrbitsAndObtainROIMetrics() {

        AbsoluteDate endDate = Utils.stamp2AD(appConfig.endDate(), DataContext.getDefault().getTimeScales().getUTC());

        List<TimedMetricsRecord> roiAccessStatistics = new ArrayList<>();

        for (AbsoluteDate date = Utils.stamp2AD(appConfig.startDate()); date.compareTo(endDate) <= 0; date = date.shiftedBy(appConfig.timeStep())) {

            long timeElapsed = TimeUtils.stamp2unix(date.toString()) - TimeUtils.stamp2unix(appConfig.startDate());

            // Obtain the starting non-euclidean FOVs
            List<AccessRegion> nonEuclideanAccessRegions = computeAccessRegionsAt(satelliteList, simulation, date, timeElapsed);

            // Check satellites within the ROI
            computeTrappedParticleEstimation(roiAccessStatistics, nonEuclideanAccessRegions, timeElapsed);

        }

        return roiAccessStatistics;
    }

    private void computeTrappedParticleEstimation(List<TimedMetricsRecord> roiAccessStatistics, List<AccessRegion> nonEuclideanAccessRegions, long timestamp) {
        TimedMetricsRecord timedMetricsRecord = new TimedMetricsRecord(timestamp, satelliteList.size());
        for (AccessRegion region : nonEuclideanAccessRegions) {
            double[] geoCoordinates = {region.getSspLat(), region.getSspLon()};
            // Get the TPO value
            double euclideanDistance = Double.MAX_VALUE;
            double tpoValue = 0;
            for (double[] coordinate : nonEuclideanROI) {
                double currentEuclideanDistance = Math.sqrt(Math.pow(coordinate[0] - geoCoordinates[0],2)
                        + Math.pow(coordinate[1] - geoCoordinates[1],2));
                if (currentEuclideanDistance < euclideanDistance) {
                    tpoValue = coordinate[2];
                    euclideanDistance = currentEuclideanDistance;
                }
                if (currentEuclideanDistance < 2) // less than 2 degrees is fine
                    break;
            }
            timedMetricsRecord.addMetric(region.getSatId(), tpoValue);
        }
        roiAccessStatistics.add(timedMetricsRecord);
    }

    private void computeSatellitesWithinROI(List<TimedMetricsRecord> roiAccessStatistics, List<AccessRegion> nonEuclideanAccessRegions, long timestamp) {

        TimedMetricsRecord timedMetricsRecord = new TimedMetricsRecord(timestamp, satelliteList.size());

        List<double[]> euclideanRoi = Transformations.toEuclideanPlane(nonEuclideanROI, 0.0, 0.0);
        for (AccessRegion region : nonEuclideanAccessRegions) {
            double[] euclideanSSP = Transformations.toStereo(Math.toRadians(region.getSspLat()), Math.toRadians(region.getSspLon()), 0.0, 0.0);
            boolean isTheSatelliteWithinTheROI = PolygonOperator.pointInPolygon(euclideanRoi, euclideanSSP);
            timedMetricsRecord.addMetric(region.getSatId(), isTheSatelliteWithinTheROI ? 1 : 0);
        }

        roiAccessStatistics.add(timedMetricsRecord);

    }

    @SuppressWarnings("unused")
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

    @SuppressWarnings("unused")
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

        // Timekeeping
        AbsoluteDate startDate = Utils.stamp2AD(appConfig.startDate());
        AbsoluteDate endDate = Utils.stamp2AD(appConfig.endDate());
        long startTimestamp = TimeUtils.stamp2unix(appConfig.startDate());

        List<AccessAreaPolygon> roiIntersections = new ArrayList<>();
        List<AccessAreaPolygon> roiUnions = new ArrayList<>();
        this.constellationCoverageTimeSeries = new ArrayList<>();
        AtomicInteger changes = new AtomicInteger();

        Log.info("Analyzing ROI coverage");
        for (AbsoluteDate t = startDate; t.compareTo(endDate) <= 0; t = t.shiftedBy(appConfig.timeStep())) {

            long timeElapsed = TimeUtils.stamp2unix(t.toString()) - startTimestamp;

            // Group regions by number of satellites on sight, for this particular time step
            Map<Integer, CopyOnWriteArrayList<AccessAreaPolygon>> byAssetsInSight = mapByNOfAssets(accessAreaPolygons, timeElapsed);

            // Guard against missing k=1 entry
            if (!byAssetsInSight.containsKey(1) || byAssetsInSight.get(1).isEmpty()) {
                Log.warn("No single-satellite AAPs at t=" + timeElapsed + "; skipping timestep.");
                constellationCoverageTimeSeries.add(new TimedMetricsRecord(timeElapsed, satelliteList.size()));
                statistics.add(timeElapsed + ",".repeat(satelliteList.size()));
                continue;
            }

            double[] surfaceValues = new double[satelliteList.size()];

            // Starting euclidean ROI
            AtomicReference<Double> roiReferenceLat = new AtomicReference<>(byAssetsInSight.get(1).get(0).getReferenceLat());
            AtomicReference<Double> roiReferenceLon = new AtomicReference<>(byAssetsInSight.get(1).get(0).getReferenceLon());
            AtomicReference<List<double[]>> euclideanROI = new AtomicReference<>(Transformations.toEuclideanPlane(nonEuclideanROI, roiReferenceLat.get(), roiReferenceLon.get()));

            TimedMetricsRecord timedMetricsRecord = new TimedMetricsRecord(timeElapsed, satelliteList.size());

            byAssetsInSight.forEach((k, aaps) -> {

                if (k < 1 || k > satelliteList.size()) {
                    Log.warn("Unexpected k=" + k + "; skipping group.");
                    return;
                }

                double referenceLat = aaps.get(0).getReferenceLat();
                double referenceLon = aaps.get(0).getReferenceLon();

                // Only if the reference changes, re-project
                if (referenceLat != roiReferenceLat.get() || referenceLon != roiReferenceLon.get()) {
                    changes.getAndIncrement();
                    roiReferenceLat.set(referenceLat);
                    roiReferenceLon.set(referenceLon);
                    euclideanROI.set(Transformations.toEuclideanPlane(nonEuclideanROI, referenceLat, referenceLon));
                }

                // Sequential stream — safe for plain ArrayList and double[]
                List<List<double[]>> unionQueue = new ArrayList<>();

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

                        List<List<double[]>> allRegions = intersection.getRegions();

                        if (allRegions.isEmpty()) {
                            return;
                        }

                        if (allRegions.size() > 1) {
                            Log.debug("Multipart intersection: " + allRegions.size() + " parts at k=" + k + ", t=" + timeElapsed);
                        }

                        for (List<double[]> eIntersectionPart : allRegions) {

                            if (eIntersectionPart.size() < 3) {
                                Log.warn("Intersection part has < 3 points; skipping.");
                                continue;
                            }

                            List<double[]> neIntersectionPart = Transformations.toNonEuclideanPlane(eIntersectionPart, referenceLat, referenceLon);

                            if (k == 1) {
                                List<Integer> gwsInSight = accessAreaPolygon.getGwsInSight();
                                if (gwsInSight != null && !gwsInSight.isEmpty()) {
                                    timedMetricsRecord.addMetric(gwsInSight.get(0),
                                            geographicTools.computeNonEuclideanSurface(neIntersectionPart));
                                }
                            }

                            unionQueue.add(eIntersectionPart);

                            AccessAreaPolygon intersectionAccessAreaPolygon = new AccessAreaPolygon(timeElapsed, k, accessAreaPolygon.getGwsInSight(), neIntersectionPart,
                                    neIntersectionPart.stream().map(pair -> pair[0]).toList(),
                                    neIntersectionPart.stream().map(pair -> pair[1]).toList());
                            roiIntersections.add(intersectionAccessAreaPolygon);
                        }

                    } catch (RuntimeException e) {
                        Log.error("Error trying to intersect the following polygon at k=" + k + ", t=" + timeElapsed + ": " + e.getMessage());
                        accessAreaPolygon.getGeoCoordinates().forEach(c -> Log.error(c[0] + "," + c[1]));
                        Log.error("### WITH ###");
                        nonEuclideanROI.forEach(c -> Log.error(c[0] + "," + c[1]));
                        Log.error(e.getLocalizedMessage());
                    }
                });

                if (!unionQueue.isEmpty()) {
                    try {
                        Polygon union = polygonOperator.polyUnion(unionQueue);
                        union = polygonOperator.polyUnion(union.getRegions()); // Second pass for clipped polygons

                        for (List<double[]> region : union.getRegions()) {
                            if (region.size() < 3) {
                                Log.warn("Union region has < 3 points at k=" + k + "; skipping.");
                                continue;
                            }
                            List<double[]> neIntersection = Transformations.toNonEuclideanPlane(region, referenceLat, referenceLon);
                            surfaceValues[k - 1] += geographicTools.computeNonEuclideanSurface(neIntersection);

                            AccessAreaPolygon unionAccessAreaPolygon = new AccessAreaPolygon(timeElapsed, k, null, neIntersection,
                                    neIntersection.stream().map(pair -> pair[0]).toList(),
                                    neIntersection.stream().map(pair -> pair[1]).toList());
                            roiUnions.add(unionAccessAreaPolygon);
                        }
                    } catch (NullPointerException e) {
                        Log.error("Regions empty?: " + unionQueue.isEmpty());
                        Log.error("Regions size?: " + unionQueue.size());
                        Log.error(e.getLocalizedMessage());
                    }
                }
            });

            StringBuilder sb = new StringBuilder(timeElapsed + "");

            for (double surface : surfaceValues) {
                sb.append(",");
                double percentage = (surface / roiSurface) * 100D;
                sb.append(percentage);
            }

            timedMetricsRecord.scale(1 / roiSurface);
            constellationCoverageTimeSeries.add(timedMetricsRecord);
            statistics.add(sb.toString());

        }
        Log.info("Ending ROI coverage analysis");

        Log.debug("Changes: " + changes);

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