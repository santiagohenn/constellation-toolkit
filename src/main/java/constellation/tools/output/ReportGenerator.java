package constellation.tools.output;

import com.google.gson.Gson;

import java.io.FileWriter;
import java.io.IOException;
import java.util.List;

public class ReportGenerator {

    private static final String JSON_EXTENSION = ".json";
    private static final String CSV_EXTENSION = ".csv";
    private final String outputPath;

    public ReportGenerator(String outputPath) {
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




}
