import constellation.tools.ConstellationCoverageComputer;
import constellation.tools.geometry.GeographicTools;
import satellite.tools.assets.entities.Satellite;

import java.util.List;

public class Main {

    public static void main(String[] args) {

        List<double[]> ROI = GeographicTools.computeSphericalCap(8, 0.0, -8.0, 150);
        ConstellationCoverageComputer constellationCoverageComputer = new ConstellationCoverageComputer("C:\\Users\\Santi\\Desktop\\Doctorado\\d3co_analysis\\config.d3co.properties");
        constellationCoverageComputer.setROI(ROI);

        int step = 0;
        double centralRAAN = 340;
        double centralAnomaly = 285;
        while (step <= 18) {
            List<Satellite> satelliteList = constellationCoverageComputer.getSatelliteList();
            satelliteList.get(0).getElements().setRightAscension(centralRAAN - step);
            satelliteList.get(0).getElements().setAnomaly(centralAnomaly + step);
            satelliteList.get(1).getElements().setAnomaly(centralAnomaly + step);
            satelliteList.get(2).getElements().setRightAscension(centralRAAN + step);
            satelliteList.get(2).getElements().setAnomaly(centralAnomaly + step);
            satelliteList.get(3).getElements().setRightAscension(centralRAAN - step);
            satelliteList.get(5).getElements().setRightAscension(centralRAAN + step);
            satelliteList.get(6).getElements().setAnomaly(centralAnomaly - step);
            satelliteList.get(6).getElements().setRightAscension(centralRAAN - step);
            satelliteList.get(7).getElements().setAnomaly(centralAnomaly - step);
            satelliteList.get(8).getElements().setRightAscension(centralRAAN + step);
            satelliteList.get(8).getElements().setAnomaly(centralAnomaly - step);
            constellationCoverageComputer.run();
            step++;
        }

    }

}
