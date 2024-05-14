package constellation.tools.utilities;

import org.apache.commons.configuration2.Configuration;
import org.apache.commons.configuration2.JSONConfiguration;
import org.apache.commons.configuration2.builder.fluent.Configurations;
import org.apache.commons.configuration2.ex.ConfigurationException;

import java.io.File;
import java.util.Locale;

public record AppConfig(
        String orekitPath,
        String startDate,
        long unixStartDate,
        String endDate,
        double timeStep,
        String outputPath,
        String satellitesFile,
        String roiPath,
        String positionsPath,
        boolean saveSnapshot,
        long snapshot,
        boolean debugMode,
        boolean propagateInternally,
        boolean saveEuclidean,
        boolean saveGeographic,
        double visibilityThreshold,
        int polygonSegments,
        double polygonEpsilon,
        double lambdaExclusion,
        int maxSubsetSize
) {
    public AppConfig(String configFilePath) throws ConfigurationException {
        this(loadConfiguration(configFilePath));
    }

    private AppConfig(Configuration config) {
        this(
                config.getString("orekit_data_path"),
                config.getString("start_date"),
                TimeUtils.stamp2unix(config.getString("start_date")),
                config.getString("end_date"),
                config.getDouble("time_step"),
                config.getString("output_path"),
                config.getString("satellites_file"),
                config.getString("roi_path"),
                config.getString("positions_path"),
                shouldSaveSnapshot(config.getString("save_snapshot")),
                shouldSaveSnapshot(config.getString("save_snapshot")) ? config.getLong("snapshot") : 0L,
                config.getBoolean("debug_mode"),
                config.getBoolean("propagate_internally"),
                config.getBoolean("save_euclidean"),
                config.getBoolean("save_geographic"),
                config.getDouble("visibility_threshold"),
                config.getInt("polygon_segments"),
                config.getDouble("polygon_epsilon"),
                config.getDouble("lambda_exclusion"),
                config.getInt("max_subset_size")
        );
    }

    private static Configuration loadConfiguration(String configFilePath) throws ConfigurationException {
        Configurations configs = new Configurations();
        if (configFilePath.toLowerCase(Locale.ROOT).endsWith(".json")) {
            return configs.fileBased(JSONConfiguration.class, new File(configFilePath));
        }
        return configs.properties(configFilePath);
    }

    private static boolean shouldSaveSnapshot(String saveSnapshot) {
        return saveSnapshot != null && !saveSnapshot.trim().isEmpty() && !saveSnapshot.equalsIgnoreCase("false");
    }

}
