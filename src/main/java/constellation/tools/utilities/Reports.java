package constellation.tools.utilities;

import com.google.gson.Gson;
import constellation.tools.geometry.OblateAccessRegion;
import satellite.tools.structures.Ephemeris;
import satellite.tools.utils.Log;

import java.io.*;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

public class Reports {

    public static final String JSON_EXTENSION = ".json";
    public static final String CSV_EXTENSION = ".csv";
    private final String outputPath;

    public Reports(String outputPath) {
        if (!outputPath.endsWith("/")) {
            outputPath = outputPath + "/";
        }
        this.outputPath = outputPath;
    }

    /**
     * Saves an Object to a JSON interpretation
     * @param obj the object to be stored
     * @param name the name of the file
     */
    public void saveAsJSON(Object obj, String name) {
        Gson gson = new Gson();
        String json = gson.toJson(obj);
        saveString(json, outputPath + name + JSON_EXTENSION);
    }

    /**
     * Saves a String List to a new file, each entry as a line
     * @param entryList the String to be stored
     * @param name the name to the file
     */
    public boolean saveAsCSV(List<String> entryList, String name) {

        try (FileWriter writer = new FileWriter(outputPath + name + CSV_EXTENSION)) {
            for (String entry : entryList) {
                writer.write(entry + "\n");
            }
        } catch (IOException e) {
            e.printStackTrace();
            return false;
        }
        return true;
    }

    /**
     * Saves a String to a file depicted in path
     * @param entry the String to be stored
     * @param path the path to the file
     */
    public boolean saveString(String entry, String path) {
        try (FileWriter writer = new FileWriter(path)) {
            writer.write(entry);
        } catch (IOException e) {
            e.printStackTrace();
            return false;
        }
        return true;
    }

    /**
     * Saves a String List to a file depicted in path
     * @param entries the String to be stored
     * @param name the path to the file
     */
    public boolean saveAs(List<String> entries, String name) {
        try (FileWriter writer = new FileWriter(outputPath + name)) {
            for (String entry : entries) {
                writer.write(entry + '\n');
            }
        } catch (IOException e) {
            e.printStackTrace();
            return false;
        }
        return true;
    }


    public List<Map<Long, Ephemeris>> positionsFromPath(String path, int nOfSatellites) {

        List<Map<Long, Ephemeris>> constellation = new ArrayList<>();

        for (int nSat = 0; nSat < nOfSatellites; nSat++) {
            Map<Long, Ephemeris> positions = new LinkedHashMap<>();
            var file = new File(path + "S" + nSat + "" + Reports.CSV_EXTENSION);
            try (var fr = new FileReader(file); var br = new BufferedReader(fr)) {
                String line;
                int id = 0;
                while ((line = br.readLine()) != null) {
                    if (!line.startsWith("//") && line.length() > 0) {
                        var data = line.split(",");
                        long time = Long.parseLong(data[0]);
                        double x = round(Double.parseDouble(data[1]));
                        double y = round(Double.parseDouble(data[2]));
                        double z = round(Double.parseDouble(data[3]));
                        Ephemeris eph = new Ephemeris(time, x, y, z);
                        double[] ssp = OblateAccessRegion.ecef2llaD(x, y, z);
                        eph.setSSP(ssp[0], ssp[1], ssp[2]);
                        positions.put(time, eph);
                    }
                }
            } catch (FileNotFoundException e) {
                Log.warn("Unable to find file: " + file);
                e.printStackTrace();
            } catch (IOException e) {
                Log.error("IOException: " + file);
                e.printStackTrace();
            }
            constellation.add(positions);
        }
        return constellation;
    }

    private double round(double num) {
        return Math.round(num * 1000000D) / 1000000D;
    }

}
