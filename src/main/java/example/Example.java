package example;

import java.util.List;

import constellation.tools.ConstellationCoverageComputer;
import constellation.tools.geometry.GeographicTools;
import constellation.tools.utilities.FileUtils;
import satellite.tools.assets.entities.Satellite;
import satellite.tools.utils.Utils;

public class Example {

    public static void main(String[] args) {
        
        ConstellationCoverageComputer constellationCoverageComputer = new ConstellationCoverageComputer("C:/Projects/lora-simulator/src/test/resources/config.d3co.properties");
        
        // List<double[]> ROI = GeographicTools.computeSphericalCap(-54.0, -58.0, 7.1946, 150);
        
        // Load from file
        List<double[]> ROI = FileUtils.file2DoubleList("C:/Projects/lora-simulator/inputs/rois/south_america.csv");
        
        GeographicTools geographicTools = new GeographicTools();
        double roiSurface = geographicTools.computeNonEuclideanSurface(ROI);
        System.out.println("ROI surface: " + roiSurface + " km^2");

        constellationCoverageComputer.setOutputPath("C:/projects/results/optimization_debug/");
        List<Satellite> satelliteList = Utils.satellitesFromFile("C:\\Projects\\lora-simulator\\inputs\\constellation_sa.csv");

        constellationCoverageComputer.setSatelliteList(satelliteList);
        constellationCoverageComputer.setROI(ROI);

        constellationCoverageComputer.run();
        System.out.println("Simulation completed. " + "Hash: " + constellationCoverageComputer.getSimulationHash());

    }
    
}
