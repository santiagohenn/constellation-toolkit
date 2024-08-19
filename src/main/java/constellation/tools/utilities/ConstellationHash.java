package constellation.tools.utilities;

import org.jline.utils.Log;
import satellite.tools.assets.entities.Satellite;
import satellite.tools.structures.OrbitalElements;

import java.io.UnsupportedEncodingException;
import java.nio.charset.StandardCharsets;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

public class ConstellationHash {

    private static String sha256(String input) {
        try {
            MessageDigest digest = MessageDigest.getInstance("SHA-256");
            byte[] hashBytes = digest.digest(input.getBytes(StandardCharsets.UTF_8));
            StringBuilder hexString = new StringBuilder();
            for (byte b : hashBytes) {
                String hex = Integer.toHexString(0xff & b);
                if (hex.length() == 1) {
                    hexString.append('0');
                }
                hexString.append(hex);
            }
            return hexString.toString();
        } catch (NoSuchAlgorithmException e) {
            Log.error("Something terrible happened, probably with your satellite's elements file");
            throw new RuntimeException(e);
        }
    }

    public static String hash(List<Satellite> constellation) {

        StringBuilder sb = new StringBuilder();
        int[] elementsSum = new int[6];

        for (Satellite satellite : constellation) {
            OrbitalElements elements = satellite.getElements();
            elementsSum[0] += truncateToThreeDecimals(elements.getSemiMajorAxis());
            elementsSum[1] += truncateToThreeDecimals(elements.getEccentricity());
            elementsSum[2] += truncateToThreeDecimals(elements.getInclination());
            elementsSum[3] += truncateToThreeDecimals(elements.getRightAscension());
            elementsSum[4] += truncateToThreeDecimals(elements.getArgOfPerigee());
            elementsSum[5] += truncateToThreeDecimals(elements.getAnomaly());
        }

        for (int parameter : elementsSum) {
            sb.append(parameter);
        }

        return sha256(sb.toString());
    }

    public static String hash2(List<Satellite> constellation) {

        StringBuilder sb = new StringBuilder();
        // int[] elementsSum = new int[6];
        List<StringBuilder> elements_concatenated = new ArrayList<>();
        for (int i = 0; i < 6; i++) {
            elements_concatenated.add(new StringBuilder());
        }
        double k = 1;
        for (Satellite satellite : constellation) {
            OrbitalElements elements = satellite.getElements();
            elements_concatenated.get(0).append((int) (elements.getSemiMajorAxis()));
            elements_concatenated.get(1).append((int) (elements.getEccentricity() * 1000));
            elements_concatenated.get(2).append((int) (elements.getInclination() * 1000));
            elements_concatenated.get(3).append((int) (elements.getRightAscension() * 1000));
            elements_concatenated.get(4).append((int) (elements.getArgOfPerigee() * 1000));
            elements_concatenated.get(5).append((int) (elements.getAnomaly() * 1000));
        }

        Collections.sort(elements_concatenated);
        elements_concatenated.forEach(sb::append);
        String stringToHash = sb.toString().strip();

        return sha256(stringToHash);
    }

    private static double truncateToThreeDecimals(double value) {
        return Math.floor(value * 1000) / 1000.0;
    }

}
