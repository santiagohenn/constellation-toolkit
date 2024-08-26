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
                nonEuclideanROI = fileUtils.file2DoubleList(appConfig.roiPath());
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
            Log.debug("Reading positions from directory: " + appConfig.positionsPath());
            constellation = fileUtils.positionsFromPath(appConfig.positionsPath(), satelliteList.size());
        }

        List<AccessAreaPolygon> nonEuclideanAccessAreaPolygons = propagateAccessAreaPolygons();
        analyzeROICoverage(nonEuclideanAccessAreaPolygons);

    }

    private List<AccessAreaPolygon> propagateAccessAreaPolygons() {

        AbsoluteDate startDate = Utils.stamp2AD(appConfig.startDate());
        AbsoluteDate endDate = Utils.stamp2AD(appConfig.endDate(), DataContext.getDefault().getTimeScales().getUTC());

        // We compute a Utility "List of Lists", containing all possible overlapping combinations between regions.
        Combination comb = new Combination(satelliteList.size(), appConfig.maxSubsetSize());
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
        if (coordinates.size() <= 2) {
            Log.debug("Less than 2 coordinates in polygon intersection");
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
        List<TimedMetricsRecord> timeSeriesData = new ArrayList<>();
        AtomicInteger changes = new AtomicInteger();

        Log.info("Analyzing ROI coverage");
        for (AbsoluteDate t = startDate; t.compareTo(endDate) <= 0; t = t.shiftedBy(appConfig.timeStep())) {

            long timeElapsed = TimeUtils.stamp2unix(t.toString()) - startTimestamp;

            // Group regions by number of satellites on sight, for this particular time step
            Map<Integer, CopyOnWriteArrayList<AccessAreaPolygon>> byAssetsInSight = mapByNOfAssets(accessAreaPolygons, timeElapsed);

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
                aaps.parallelStream().forEach(accessAreaPolygon -> {

                    List<double[]> eIntersection = new ArrayList<>();

                    try {
                        List<List<double[]>> polygonAndROI = Arrays.asList(euclideanROI.get(),
                                Transformations.toEuclideanPlane(accessAreaPolygon.getGeoCoordinates(),
                                        referenceLat, referenceLon));
                        Polygon intersection = polygonOperator.polyIntersect(polygonAndROI);
                        eIntersection = intersection.getRegions().isEmpty() ? eIntersection : intersection.getRegions().get(0);

                        if (intersection.getRegions().size() > 1) {
                            Log.warn("More than 1 intersection!!!");
                        }

                    } catch (RuntimeException e) {
                        Log.error("Error trying to intersect the following polygon: ");
                        accessAreaPolygon.getGeoCoordinates().forEach(c -> Log.error(c[0] + "," + c[1]));
                        Log.error("### WITH ###");
                        nonEuclideanROI.forEach(c -> Log.error(c[0] + "," + c[1]));
                        e.printStackTrace();
                    }

                    if (eIntersection.size() >= 3) {
                        unionQueue.add(eIntersection);
                    } else if (!eIntersection.isEmpty()) {
                        Log.warn("Intersection with less than 3 points");
                    }

                    List<double[]> neIntersection = Transformations.toNonEuclideanPlane(eIntersection, referenceLat, referenceLon);

                    // Surface for K = 1 (intersection seen by 1 GW):
                    if (k == 1) {
                        timedMetricsRecord.addMetric(accessAreaPolygon.getGwsInSight().get(0), geographicTools.computeNonEuclideanSurface(neIntersection));
                    }

                    AccessAreaPolygon intersectionAccessAreaPolygon = new AccessAreaPolygon(timeElapsed, k, accessAreaPolygon.getGwsInSight(), neIntersection,
                            neIntersection.stream().map(pair -> pair[0]).toList(),
                            neIntersection.stream().map(pair -> pair[1]).toList());
                    roiIntersections.add(intersectionAccessAreaPolygon);

                });

                if (!unionQueue.isEmpty()) {
                    try {
                        Polygon union = polygonOperator.polyUnion(unionQueue);
                        union = polygonOperator.polyUnion(union.getRegions());  // Second union if some polygons got clipped out
                        union.getRegions().forEach(region -> {
                            List<double[]> neIntersection = Transformations.toNonEuclideanPlane(region, referenceLat, referenceLon);
                            surfaceValues[k - 1] = surfaceValues[k - 1] + geographicTools.computeNonEuclideanSurface(neIntersection);
                            AccessAreaPolygon unionAccessAreaPolygon = new AccessAreaPolygon(timeElapsed, k, null, neIntersection,
                                    neIntersection.stream().map(pair -> pair[0]).toList(),
                                    neIntersection.stream().map(pair -> pair[1]).toList());
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
                double percentage = (surface / roiSurface) * 100D;
                sb.append(percentage);
            }

            timedMetricsRecord.scale(1 / roiSurface);
            timeSeriesData.add(timedMetricsRecord);
            statistics.add(sb.toString());

        }
        Log.info("Ending ROI coverage analysis");

        Log.debug("Changes: " + changes);

        if (appConfig.saveGeographic() && appConfig.saveSnapshot()) {
            saveAAPsAt(roiIntersections, "snapshot_aaps_intersection_" + simulationHash, appConfig.snapshot());
            saveAAPsAt(roiUnions, "snapshot_aaps_union_" + simulationHash, appConfig.snapshot());
        }

        fileUtils.saveAsCSV(statistics, "coverage_" + simulationHash);
        fileUtils.saveAsJSON(timeSeriesData, "surface_metrics_" + simulationHash);

    }

    private Map<Integer, CopyOnWriteArrayList<AccessAreaPolygon>> mapByNOfAssets(List<AccessAreaPolygon> accessAreaPolygons, long timeElapsed) {

        Map<Integer, CopyOnWriteArrayList<AccessAreaPolygon>> byAssetsInSight = new LinkedHashMap<>(appConfig.maxSubsetSize());
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

        for (int nOfGw = 1; nOfGw <= appConfig.maxSubsetSize(); nOfGw++) {
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
                e.printStackTrace();
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

}