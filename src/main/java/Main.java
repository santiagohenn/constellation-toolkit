import constellation.tools.ConstellationCoverageComputer;
import constellation.tools.geometry.GeographicTools;

import java.util.List;

public class Main {

    public static void main(String[] args) {

        List<double[]> ROI = GeographicTools.computeSphericalCap(7.1946, -21.0, -58.0, 200);
        ConstellationCoverageComputer constellationCoverageComputer = new ConstellationCoverageComputer("C:\\Users\\Santi\\Desktop\\Doctorado\\D3CO analysis\\config.properties");
        constellationCoverageComputer.setROI(ROI);
        constellationCoverageComputer.run();

    }

}
