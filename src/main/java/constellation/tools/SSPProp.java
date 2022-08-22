package constellation.tools;

import constellation.tools.geometry.ConstellationSSPs;
import constellation.tools.geometry.FOV;
import constellation.tools.math.Pair;
import constellation.tools.reports.ReportGenerator;
import org.orekit.time.AbsoluteDate;
import satellite.tools.Simulation;
import satellite.tools.assets.entities.Satellite;
import satellite.tools.structures.Ephemeris;
import satellite.tools.utils.Utils;

import java.util.ArrayList;
import java.util.List;
import java.util.Properties;

public class SSPProp {

    private static final Properties prop = Utils.loadProperties("config.properties");
    private static final String orekitPath = (String) prop.get("orekit_data_path");
    private final String START_DATE = (String) prop.get("start_date");
    private final String END_DATE = (String) prop.get("end_date");
    private final double TIME_STEP = Double.parseDouble((String) prop.get("time_step"));
    private final String OUTPUT_PATH = (String) prop.get("output_path");
    private final String SATELLITES_FILE = (String) prop.get("satellites_file");
    private final List<Satellite> satelliteList = Utils.satellitesFromFile(SATELLITES_FILE);
    private final List<ConstellationSSPs> constellationSSPs = new ArrayList<>();
    ReportGenerator reportGenerator = new ReportGenerator(OUTPUT_PATH);

    public SSPProp() {

    }

    public void runSSPs() {

        Simulation simulation = new Simulation(orekitPath);

        AbsoluteDate endDate = Utils.stamp2AD(END_DATE);
        AbsoluteDate startDate = Utils.stamp2AD(START_DATE);
        AbsoluteDate pointerDate = startDate;
        double scenarioDuration = endDate.durationFrom(startDate);

        while (pointerDate.compareTo(endDate) <= 0) {

            long timeSinceStart = Utils.stamp2unix(pointerDate.toString()) - Utils.stamp2unix(START_DATE);
            updateProgressBar(pointerDate.durationFrom(startDate), scenarioDuration);

            // Obtain the starting non-euclidean FOVs and their surface value
            List<FOV> nonEuclideanFOVs = computeSSPsAt(satelliteList, simulation, pointerDate);

            // Store the SSPs (per Guido's request)
            List<Pair> SSPs = new ArrayList<>();

            nonEuclideanFOVs.forEach(FOV -> SSPs.add(new Pair(FOV.getSspLat(), FOV.getSspLon())));
            constellationSSPs.add(new ConstellationSSPs(timeSinceStart, SSPs));

            // Advance to the next time step
            pointerDate = pointerDate.shiftedBy(TIME_STEP);

        }

        reportGenerator.saveAsJSON(constellationSSPs, "ConstellationSSPs");

    }

    /**
     * This method takes the satellite list, propagates orbits to the specified date and computes the corresponding
     * access area or FOV polygon as a Region object for each one.
     *
     * @param satelliteList a list of Satellite objects
     * @param date          an AbsoluteDate object
     * @return a List of Regions
     **/
    private List<FOV> computeSSPsAt(List<Satellite> satelliteList, Simulation simulation, AbsoluteDate date) {

        List<FOV> FOVList = new ArrayList<>();

        for (Satellite satellite : satelliteList) {
            simulation.setSatellite(satellite);
            Ephemeris ephemeris = simulation.computeEphemerisKm(date);
            FOV FOV = new FOV(satellite.getId(), ephemeris.getLatitude(), ephemeris.getLongitude(), null);
            FOVList.add(FOV);
        }

        return FOVList;

    }

    /**
     * Rudimentary progress bar.
     * **/
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

}
