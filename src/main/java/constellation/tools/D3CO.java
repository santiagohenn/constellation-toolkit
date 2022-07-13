package constellation.tools;

import com.menecats.polybool.Epsilon;
import com.menecats.polybool.PolyBool;
import com.menecats.polybool.models.Polygon;
import constellation.tools.geometry.AAP;
import constellation.tools.geometry.FOV;
import constellation.tools.geometry.Geo;
import constellation.tools.math.Combination;
import constellation.tools.math.Transformations;
import constellation.tools.output.ReportGenerator;
import org.orekit.time.AbsoluteDate;
import satellite.tools.Simulation;
import satellite.tools.assets.entities.Satellite;
import satellite.tools.structures.Ephemeris;
import satellite.tools.utils.Log;
import satellite.tools.utils.Utils;

import java.util.*;
import java.util.stream.Collectors;

import static com.menecats.polybool.helpers.PolyBoolHelper.epsilon;

/**
 * Dynamic Constellation Coverage Computer (D3CO)
 **/
public class D3CO {

    private static final Properties prop = Utils.loadProperties("config.properties");
    private static final String orekitPath = (String) prop.get("orekit_data_path");
    private final String START_DATE = (String) prop.get("start_date");
    private final long UNIX_START_DATE = Utils.stamp2unix(START_DATE);
    private final String END_DATE = (String) prop.get("end_date");
    private final double TIME_STEP = Double.parseDouble((String) prop.get("time_step"));
    private final String OUTPUT_PATH = (String) prop.get("output_path");
    private final String SATELLITES_FILE = (String) prop.get("satellites_file");
    private final String ROI_PATH = (String) prop.get("roi_path");
    private final boolean DEBUG = Boolean.parseBoolean((String) prop.get("debug_mode"));
    private final long SNAPSHOT = Long.parseLong((String) prop.get("snapshot"));
    private final double VISIBILITY_THRESHOLD = Double.parseDouble((String) prop.get("visibility_threshold"));
    private final double POLYGON_SEGMENTS = Double.parseDouble((String) prop.get("polygon_segments"));
    private final int MAX_SUBSET_SIZE = Integer.parseInt((String) prop.get("max_subset_size"));

    private final List<Satellite> satelliteList = Utils.satellitesFromFile(SATELLITES_FILE);
    private final List<String> statistics = new ArrayList<>();

    ReportGenerator reportGenerator = new ReportGenerator(OUTPUT_PATH);

    /**
     * Default constructor
     **/
    public D3CO() {


    }

    public void run() {

        Simulation simulation = new Simulation(orekitPath);

        AbsoluteDate endDate = Utils.stamp2AD(END_DATE);
        AbsoluteDate startDate = Utils.stamp2AD(START_DATE);
        AbsoluteDate pointerDate = startDate;
        double scenarioDuration = endDate.durationFrom(startDate);

        // We compute a Utility "List of Lists", containing all possible overlapping combinations between regions.
        Combination comb = new Combination(satelliteList.size(), MAX_SUBSET_SIZE);
        final List<List<Integer>> combinationsList = comb.computeCombinations();

        List<AAP> nonEuclideanAAPs = new ArrayList<>();
        List<AAP> euclideanAAPs = new ArrayList<>();

        double lambdaMax = Geo.getLambdaMax(satelliteList.get(0).getElement("a"), VISIBILITY_THRESHOLD); // FIXME do I use this?

        while (pointerDate.compareTo(endDate) <= 0) {

            long timeSinceStart = Utils.stamp2unix(pointerDate.toString()) - Utils.stamp2unix(START_DATE);
            updateProgressBar(pointerDate.durationFrom(startDate), scenarioDuration);

            // Obtain the starting non-euclidean FOVs and their surface value
            List<FOV> nonEuclideanFOVs = computeFOVsAt(satelliteList, simulation, pointerDate);

            // double surfaceInKm = 0D;

            for (List<Integer> combination : combinationsList) {

                List<double[]> nonEuclideanCoordinates = new ArrayList<>();
                List<double[]> euclideanCoordinates = new ArrayList<>();

                // assemble a list of the FOVs to be intersected at this time step:
                List<FOV> FOVsToIntersect = new ArrayList<>();
                combination.forEach(regionIndex -> FOVsToIntersect.add(nonEuclideanFOVs.get(regionIndex)));

                // Get reference point for the projection
                int poleProximity = checkPoleInclusion(FOVsToIntersect, lambdaMax);
                double referenceLat = poleProximity * 90; // FOVsToIntersect.get(0).getReferenceLat();
                double referenceLon = 0; // FOVsToIntersect.get(0).getReferenceLon();

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

                euclideanAAPs.add(new AAP(timeSinceStart, combination.size(), combination, euclideanCoordinates,
                        euclideanCoordinates.stream().map(pair -> pair[0]).collect(Collectors.toList()),
                        euclideanCoordinates.stream().map(pair -> pair[1]).collect(Collectors.toList())));

            }

            // Advance to the next time step
            pointerDate = pointerDate.shiftedBy(TIME_STEP);

        }

        saveAAPs(nonEuclideanAAPs, "NEPolygons");
        saveAAPs(euclideanAAPs, "EPolygons");

        saveAAPsAt(nonEuclideanAAPs, "NEPolygons_debug", SNAPSHOT);
        saveAAPsAt(euclideanAAPs, "EPolygons_debug", SNAPSHOT);

        analyzeSurfaceCoverage(nonEuclideanAAPs);
        analyzeROICoverage(nonEuclideanAAPs);

    }

    public void analyzeROICoverage(List<AAP> AAPs) {

        statistics.clear();

        // Load ROI Data:
        List<double[]> nonEuclideanROI = Geo.file2DoubleList(ROI_PATH);
        double roiSurface = Geo.computeNonEuclideanSurface2(nonEuclideanROI);
        Log.debug("Area of ROI in Km2: " + roiSurface / 1e6);

        // TODO: Generalize for any ROI
        double referenceLat = -90;
        double referenceLon = 0;

        // FIXME: use euclidean interception! I'm just testing this!
        List<double[]> euclideanROI = Transformations.toEuclideanPlane(nonEuclideanROI, referenceLat, referenceLon); //nonEuclideanROI; //

        // Timekeeping
        AbsoluteDate startDate = Utils.stamp2AD(START_DATE);
        AbsoluteDate pointerDate = Utils.stamp2AD(START_DATE);
        AbsoluteDate endDate = Utils.stamp2AD(END_DATE);
        double scenarioDuration = endDate.durationFrom(startDate);

        // Accumulated areas by number of satellites in visibility is stored in this array (idx = number of sats, value = area) // FIXME remove eventually
        Map<Integer, Double> accumCoverage = new HashMap<>(MAX_SUBSET_SIZE);

        List<AAP> roiIntersections2 = new ArrayList<>();
        List<AAP> roiIntersections = new ArrayList<>();

        while (pointerDate.compareTo(endDate) <= 0) {

            updateProgressBar(pointerDate.durationFrom(startDate), scenarioDuration);
            accumCoverage.clear();

            long timeElapsed = Utils.stamp2unix(pointerDate.toString()) - Utils.stamp2unix(START_DATE);

            // Group regions by number of satellites on sight, for this particular timestep
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

            double[] surfaceValues = new double[satelliteList.size()];
            double[] surfaceValuesU = new double[satelliteList.size()];

            // Perform intersection of AAPs with the ROI and surface area values calculation
            // For each number of assets
            byAssetsInSight.forEach((key, value) -> {

                List<List<double[]>> toBeOr = new ArrayList<>();

                // For each AAP with this number of assets in sight
                value.forEach(aap -> {
                    // Intersections with ROI
                    List<double[]> eIntersection = intersectAndGetPolygon(euclideanROI,
                            Transformations.toEuclideanPlane(aap.getNonEuclideanCoordinates(),
                                    referenceLat, referenceLon));

                    toBeOr.add(eIntersection);

                    List<double[]> neIntersection = Transformations.toNonEuclideanPlane(eIntersection, referenceLat, referenceLon);

                    surfaceValues[key - 1] = surfaceValues[key - 1] + Geo.computeNonEuclideanSurface2(neIntersection);

                    AAP intersectionAAP = new AAP(timeElapsed, key, aap.getGwsInSight(), neIntersection,
                            neIntersection.stream().map(pair -> pair[0]).collect(Collectors.toList()),
                            neIntersection.stream().map(pair -> pair[1]).collect(Collectors.toList()));
                    roiIntersections.add(intersectionAAP);

                });

                // Union of all ROI intersections
                Polygon polygon = new Polygon(toBeOr);
                Epsilon eps = epsilon();

                try {
                    Polygon union = PolyBool.union(eps, new Polygon(), polygon);

                    union.getRegions().forEach(region -> {

                        List<double[]> neIntersection = Transformations.toNonEuclideanPlane(region, referenceLat, referenceLon);

                        surfaceValuesU[key - 1] = surfaceValuesU[key - 1] + Geo.computeNonEuclideanSurface2(neIntersection);

                        AAP unionAAP = new AAP(timeElapsed, key, null, neIntersection,
                                neIntersection.stream().map(pair -> pair[0]).collect(Collectors.toList()),
                                neIntersection.stream().map(pair -> pair[1]).collect(Collectors.toList()));
                        roiIntersections2.add(unionAAP);

                    });
                } catch (IndexOutOfBoundsException e) {
                    Log.error(e.getMessage());
                    Log.error("toBeOr size: " + toBeOr.size());
                }

            });

//            double[] percentageValues = new double[satelliteList.size()];

            // remove intersections surfaces
            for (int i = 0; i < surfaceValues.length - 1; i++) {
                for (int j = i + 1; j < surfaceValues.length; j++) {
                    surfaceValues[i] = surfaceValues[i] - surfaceValues[j];
                }
            }

            // TODO REMOVE DEBUG
            Log.info("----------------------------------");
            Log.info("Surface values with intersection: ");
            for (double surface : surfaceValues) {
                Log.info("Normalized surface [km2]: " + surface + " => " + (surface / roiSurface) * 100.00000 + "% Coverage");
            }
            Log.info("----------------------------------");
            Log.info("Surface values with union: ");
            for (double surface : surfaceValuesU) {
                Log.info("Normalized surface [km2]: " + surface + " => " + (surface / roiSurface) * 100.00000 + "% Coverage");
            }
            Log.info("----------------------------------");

//            // Transform to euclidean plane and calculate intersection of AAPs with the ROI
//            byAssetsInSight.forEach((key, value) -> {
//                List<List<double[]>> toBeOR = new ArrayList<>();
//                // intersect for each number of sats
//                value.forEach(aap -> {
//                    List<double[]> intersection = intersectAndGetPolygon(euclideanROI,
//                            Transformations.toEuclideanPlane(aap.getNonEuclideanCoordinates(),
//                                    referenceLat, referenceLon));
//
//                    if (intersection.size() >= 3) {
//                        toBeOR.add(intersection);
//                    }
//
//                    // FIXME Might need to move this up
//                    AAP intersectionAAP = new AAP(timeElapsed, key, aap.getGwsInSight(), intersection,
//                            intersection.stream().map(pair -> pair[0]).collect(Collectors.toList()),
//                            intersection.stream().map(pair -> pair[1]).collect(Collectors.toList()));
//                    roiIntersections.add(intersectionAAP);
//                    Log.debug("intersection for : " + key + " sats: " + intersection.size() + " at " + timeElapsed + " roi intersection size: " + roiIntersections.size());
//
////                    if (roiIntersections.containsKey(key)) {
////                        roiIntersections.get(key).add(intersectionAAP);
////                    } else {
////                        List<AAP> aapList = new ArrayList<>();
////                        aapList.add(intersectionAAP);
////                        roiIntersections.put(key, aapList);
////                    }
//                });
//
//                Log.debug("To be or size: " + toBeOR.size());
//                List<List<double[]>> union = new ArrayList<>();
//
//                try {
//                    union = polyUnion(euclideanROI, toBeOR);
//                } catch (IndexOutOfBoundsException e) {
//                    Log.error("----ERROR----");
//                    Log.error(e.getMessage());
//                    toBeOR.forEach(poly -> Log.error("size: " + poly.size()));
//                    Log.error("-------------");
//                }
//
//                double areaCovered = 0;
//                for (List<double[]> poly : union) {
//                    areaCovered = areaCovered + Geo.computeNonEuclideanSurface2(Transformations.toNonEuclideanPlane(poly, referenceLat, referenceLon));
//
//                    // TODO: REMOVE THIS DEBUG
//                    AAP orAAP = new AAP(timeElapsed, key, null, poly,
//                            poly.stream().map(pair -> pair[0]).collect(Collectors.toList()),
//                            poly.stream().map(pair -> pair[1]).collect(Collectors.toList()));
//                    roiIntersections2.add(orAAP);
//
//                }
//                Log.debug("surface union: " + (areaCovered));
//                Log.debug("surface roi: " + (roiSurface));
//                Log.info("coverage: " + (areaCovered / roiSurface) * 100.0 + " %");
//
//            });


//                List<double[]> roiIntersection = intersectAndGetPolygon(ROI, AAP.getNonEuclideanCoordinates());
//
//                double coverage;
//
//                if (roiIntersection.isEmpty()) {
//                    coverage = 0;
//                } else {
//                    coverage = Geo.computeNonEuclideanSurface2(roiIntersection);
//
//                    int nAssets = AAP.getnOfGwsInSight();
//                    // accumulatedAreas.putIfAbsent(nAssets, surfaceInKm2);
//                    if (accumCoverage.containsKey(nAssets)) {
//                        accumCoverage.put(nAssets, accumCoverage.get(nAssets) + coverage);
//                    } else {
//                        accumCoverage.put(nAssets, coverage);
//                    }
//                    AAP.setSurfaceInKm2(coverage);
//
//                    // Save AAPs
//                    roiIntersections2.add(new AAP(timeElapsed, AAP.getnOfGwsInSight(), AAP.getGwsInSight(), roiIntersection,
//                            roiIntersection.stream().map(pair -> pair[0]).collect(Collectors.toList()),
//                            roiIntersection.stream().map(pair -> pair[1]).collect(Collectors.toList())));
//                }
//
//            });

            // Transform accumulated coverage into percentage
//            accumCoverage.keySet().forEach(nOfAssets -> {
//                accumCoverage.put(nOfAssets, accumCoverage.get(nOfAssets) / roiSurface);
//            });
//
//            // Save the results
//            statistics.add(stringifyResults(pointerDate, accumCoverage));

            // Advance to the next time step
            pointerDate = pointerDate.shiftedBy(TIME_STEP);

        }
        saveAAPsAt(roiIntersections, "RoiDebug", SNAPSHOT);
//        saveAAPsAt(roiIntersections2, "RoiDebug2", SNAPSHOT);

//        reportGenerator.saveAsCSV(statistics, "PercentageCoverage");

    }

    public void analyzeSurfaceCoverage(List<AAP> AAPs) {

        statistics.clear();

        AbsoluteDate endDate = Utils.stamp2AD(END_DATE);
        AbsoluteDate pointerDate = Utils.stamp2AD(START_DATE);

        // Accumulated areas by number of satellites in visibility is stored in this array (idx = number of sats, value = area) // FIXME remove eventually
        Map<Integer, Double> accumulatedAreas = new HashMap<>(MAX_SUBSET_SIZE);

        while (pointerDate.compareTo(endDate) <= 0) {   // TODO while/filter can be improved

            accumulatedAreas.clear();

            long timeSinceStart = Utils.stamp2unix(pointerDate.toString()) - Utils.stamp2unix(START_DATE);
            AAPs.stream().filter(AAP -> AAP.getDate() == timeSinceStart).forEach(AAP -> {

                double surfaceInKm2 = Geo.computeNonEuclideanSurface2(AAP.getNonEuclideanCoordinates()) * 1E-6; // AAP.getSurfaceInKm2(); //
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

    // FIXME REPLACE WITH POST-ANALYSIS
    private String stringifyResults(AbsoluteDate pointerDate, Map<Integer, Double> accumulatedAreas) {

        StringBuilder sb = new StringBuilder();
        sb.append(Utils.stamp2unix(pointerDate.toString()) - UNIX_START_DATE);

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
            int proximity = Geo.checkPoleInclusion(FOV, lambdaMax);
            if (proximity != 0) {
                return proximity;
            }
        }
        return 1;


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

    // TODO: this can be improved inheriting the library's intersection capabilities, for now, we dont trust them

    /**
     * This method takes two polygons, and returns their union using the Martinez-Rueda Algorithm.
     *
     * @see <a href="https://github.com/Menecats/polybool-java">Menecats-Polybool</a>
     * @see <a href="https://www.sciencedirect.com/science/article/pii/S0965997813000379">Martinez-Rueda clipping algorithm</a>
     **/
    private List<double[]> uniteAndGetPolygon(List<double[]> polygonA, List<double[]> polygonB) {

        List<List<double[]>> regions1 = new ArrayList<>();
        regions1.add(polygonA);

        List<List<double[]>> regions2 = new ArrayList<>();
        regions2.add(polygonB);

        Polygon polyA = new Polygon(regions1);
        Polygon polyB = new Polygon(regions2);
        Polygon union = new Polygon();

        if (polyA.getRegions().get(0).size() >= 3 && polyB.getRegions().get(0).size() >= 3) {
            Epsilon eps = epsilon();
            union = PolyBool.union(eps, polyA, polyB);
        }

        if (union.getRegions().size() > 0) {
            return union.getRegions().get(0);
        } else {
            return new ArrayList<>();
        }

    }

    /**
     * This method takes a polygon (A) and a polygon list and returns the unions between A and each polygon in the list.
     *
     * @see <a href="https://github.com/Menecats/polybool-java">Menecats-Polybool</a>
     * @see <a href="https://www.sciencedirect.com/science/article/pii/S0965997813000379">Martinez-Rueda clipping algorithm</a>
     **/
    private List<List<double[]>> polyUnion(List<double[]> polygonA, List<List<double[]>> polygonList) {

        List<List<double[]>> unionsList = new ArrayList<>();

//        List<double[]> intersectedPolygon = new ArrayList<>(polygonsToIntersect.get(0));
//
//        // Obtain access polygon
//        for (List<double[]> polygon : polygonsToIntersect) {
//            if (polygonsToIntersect.indexOf(polygon) == 0) continue;
//            intersectedPolygon = intersectAndGetPolygon(intersectedPolygon, polygon);
//        }

        // ROI
        List<List<double[]>> roi = new ArrayList<>();
        List<List<double[]>> aap = new ArrayList<>();
        roi.add(polygonA);
        Polygon polyA = new Polygon(roi);

        Epsilon eps = epsilon();
        Polygon union;

        for (List<double[]> polygon : polygonList) {

            aap.clear();
            aap.add(polygon);
            Polygon polyB = new Polygon(aap);

            try {
                union = PolyBool.union(eps, polyA, polyB);

                if (union.getRegions().size() > 0) {
                    unionsList.add(union.getRegions().get(0));
                }

            } catch (RuntimeException e) {
                Log.error(e.getMessage());
                polygonList.forEach(poly -> {
                    Log.error("POLYGON SIZE: " + poly.size());
//                    poly.forEach(pair -> {
//                        Log.error(pair[0] + "," + pair[1]);
//                    });
                });
            }

        }


        return unionsList;

    }

    /**
     * This method takes a polygon A and a list of polygons L, and returns a list of intersections between A and every
     * member of L, using the Martinez-Rueda Algorithm.
     *
     * @see <a href="https://github.com/Menecats/polybool-java">Menecats-Polybool</a>
     * @see <a href="https://www.sciencedirect.com/science/article/pii/S0965997813000379">Martinez-Rueda clipping algorithm</a>
     **/
    private List<List<double[]>> intersectWithRegions(List<double[]> polygonA, List<List<double[]>> polygonList) {

        List<List<double[]>> regions1 = new ArrayList<>();
        regions1.add(polygonA);

        Polygon polyA = new Polygon(regions1);
        Polygon polyB = new Polygon(polygonList);
        Polygon intersection = new Polygon();

        if (polyA.getRegions().get(0).size() >= 3 && polygonList.size() > 0) {
            Epsilon eps = epsilon();
            intersection = PolyBool.intersect(eps, polyA, polyB);
            return intersection.getRegions();
        } else {
            return new ArrayList<>();
        }

    }


    /**
     * Over the top progress bar mainly for debugging.
     **/
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

        if ((int) progress >= 100) {
            System.out.println();
        }

    }

}