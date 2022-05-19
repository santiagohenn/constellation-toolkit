package constellation.tools;

import org.junit.Test;
import satellite.tools.Simulation;
import satellite.tools.utils.Log;

import java.util.Properties;

public class GaapTest {

    @Test
    public void SimulationTest() {

        Simulation simulation = new Simulation();
        Properties p = simulation.getProperties();
        Log.debug("Orekit data path: " + p.get("orekit_data_path"));
        Log.debug("Start date: " + p.get("start_date"));
        Log.debug("End date: " + p.get("end_date"));
        Log.debug("Time step: " + p.get("time_step"));
        Log.debug("Vis. TH: " + p.get("visibility_threshold"));
        Log.debug("th_detection: " + p.get("th_detection"));

    }

}
