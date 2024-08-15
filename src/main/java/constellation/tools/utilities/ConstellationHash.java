package constellation.tools.utilities;

import satellite.tools.assets.entities.Satellite;
import satellite.tools.structures.OrbitalElements;

import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.util.List;

public class ConstellationHash {

    private static String sha256(String input) {
        try {
            MessageDigest digest = MessageDigest.getInstance("SHA-256");
            byte[] hashBytes = digest.digest(input.getBytes());
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

    private static double truncateToThreeDecimals(double value) {
        return Math.floor(value * 1000) / 1000.0;
    }

}
