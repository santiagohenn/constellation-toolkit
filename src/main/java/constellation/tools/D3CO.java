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

import java.io.File;
import java.util.*;
import java.util.stream.Collectors;

import static com.menecats.polybool.helpers.PolyBoolHelper.epsilon;

/**
 * Dynamic Constellation Coverage Computer (D3CO)
 **/
public class D3CO {

    private static final Properties prop = Utils.loadProperties("config.properties");
    private static final String orekitPath = (String) prop.get("orekit_data_path");
    private static final File orekitFile = Utils.loadFile(orekitPath);

    private static final DataProvidersManager manager = DataContext.getDefault().getDataProvidersManager();

    static {
        manager.addProvider(new DirectoryCrawler(orekitFile));
        Log.debug("Manager loaded");
    }

    private final String START_DATE = (String) prop.get("start_date");
    private final String END_DATE = (String) prop.get("end_date");
    private final double TIME_STEP = Double.parseDouble((String) prop.get("time_step"));
    private final String OUTPUT_PATH = (String) prop.get("output_path");
    private final String SATELLITES_FILE = (String) prop.get("satellites_file");
    private final boolean DEBUG = Boolean.parseBoolean((String) prop.get("debug_mode"));
    private final long SNAPSHOT = Long.parseLong((String) prop.get("snapshot"));
    private final double VISIBILITY_THRESHOLD = Double.parseDouble((String) prop.get("visibility_threshold"));
    private final double POLYGON_SEGMENTS = Double.parseDouble((String) prop.get("polygon_segments"));
    private final int MAX_SUBSET_SIZE = Integer.parseInt((String) prop.get("max_subset_size"));

    private final List<Satellite> satelliteList = Utils.satellitesFromFile(SATELLITES_FILE);
    private final List<String> statistics = new ArrayList<>();
    private final List<ConstellationSSPs> constellationSSPs = new ArrayList<>();

    ReportGenerator reportGenerator = new ReportGenerator(OUTPUT_PATH);

    /**
     * Default constructor
     **/
    public D3CO() {


    }

    public void runSSPs() {

        AbsoluteDate endDate = Utils.stamp2AD(END_DATE);
        AbsoluteDate startDate = Utils.stamp2AD(START_DATE);
        AbsoluteDate pointerDate = startDate;
        double scenarioDuration = endDate.durationFrom(startDate);

        while (pointerDate.compareTo(endDate) <= 0) {

            long timeSinceStart = Utils.stamp2unix(pointerDate.toString()) - Utils.stamp2unix(START_DATE);
            updateProgressBar(pointerDate.durationFrom(startDate), scenarioDuration);

            // Obtain the starting non-euclidean FOVs and their surface value
            List<FOV> nonEuclideanFOVs = computeFOVsAt(satelliteList, pointerDate);

            // Store the SSPs (per Guido's request)
            List<Pair> SSPs = new ArrayList<>();

            nonEuclideanFOVs.forEach(FOV -> SSPs.add(new Pair(FOV.getReferenceLat(), FOV.getReferenceLon())));
            constellationSSPs.add(new ConstellationSSPs(timeSinceStart, SSPs));

            // Advance to the next time step
            pointerDate = pointerDate.shiftedBy(TIME_STEP);

        }

        saveSSPs(constellationSSPs);

        reportGenerator.saveAsCSV(statistics, "stats");

    }

    public void run() {

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
            List<FOV> nonEuclideanFOVs = computeFOVsAt(satelliteList, pointerDate);

            // Accumulated areas by number of satellites in visibility is stored in this array (idx = number of sats, value = area) // FIXME remove eventually
            Map<Integer, Double> accumulatedAreas = new HashMap<>(MAX_SUBSET_SIZE);

            double surfaceInKm = 0D;

            for (List<Integer> combination : combinationsList) {

                List<Pair> nonEuclideanCoordinates = new ArrayList<>();
                List<Pair> euclideanCoordinates = new ArrayList<>();

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
                    surfaceInKm = neFov.getSurfaceKm2();
                    nonEuclideanCoordinates = Transformations.doubleList2pairList(nonEuclideanFOVs.get(fovIdx).getPolygonCoordinates());
                    euclideanCoordinates = Transformations.doubleList2pairList(Transformations
                            .toEuclideanPlane(neFov.getPolygonCoordinates(), referenceLat, referenceLon));

                } else if (checkDistances(combination, nonEuclideanFOVs, lambdaMax)) {    // FIXME optimize transforming every starting AAP to stereographic


                    List<List<double[]>> polygonsToIntersect = new ArrayList<>();
                    FOVsToIntersect.forEach(FOV -> polygonsToIntersect.add(Transformations.toEuclideanPlane(FOV.getPolygonCoordinates(), referenceLat, referenceLon)));

                    List<double[]> resultingPolygon = new ArrayList<>(polygonsToIntersect.get(0));

                    for (List<double[]> polygon : polygonsToIntersect) {
                        if (polygonsToIntersect.indexOf(polygon) == 0) continue;

                        resultingPolygon = intersectAndGetPolygon(resultingPolygon, polygon);

                    }
                    euclideanCoordinates = Transformations.doubleList2pairList(resultingPolygon);

                    List<double[]> nonEuclideanIntersection = Transformations.toNonEuclideanPlane(resultingPolygon, referenceLat, referenceLon);
                    nonEuclideanCoordinates = Transformations.doubleList2pairList(nonEuclideanIntersection);
                    Log.debug("nonEuclideanIntersection: " + nonEuclideanIntersection.size() + " - nonEuclideanSize: " + nonEuclideanCoordinates.size());

                    if (timeSinceStart == SNAPSHOT) {
                        nonEuclideanCoordinates.forEach(System.out::println);
                    }

                    surfaceInKm = Geo.computeNonEuclideanSurface(nonEuclideanCoordinates) * 1E-6;

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

        saveAAPs(nonEuclideanAAPs, "NEPolygons");
        saveAAPs(euclideanAAPs, "EPolygons");

        saveAAPsAt(nonEuclideanAAPs, "NEPolygons_debug", SNAPSHOT);
        saveAAPsAt(euclideanAAPs, "EPolygons_debug", SNAPSHOT);

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
    private List<FOV> computeFOVsAt(List<Satellite> satelliteList, AbsoluteDate date) {

        Simulation simulation = new Simulation();
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
        return 0;
    }

    /**
     * ... Polybool intersection from list of double[]
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

        return intersection.getRegions().get(0);

    }

}