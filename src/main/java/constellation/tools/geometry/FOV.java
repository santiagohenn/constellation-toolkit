package constellation.tools.geometry;

import java.util.List;

public class FOV {

    private final double referenceLat;
    private final double referenceLon;
    private final int satId;
    private List<double[]> polygonCoordinates;
    private double surface;

    public FOV(int satId, double centerLat, double centerLong, List<double[]> polygon) {
        this.satId = satId;
        this.referenceLat = centerLat;
        this.referenceLon = centerLong;
        this.polygonCoordinates = polygon;
    }

    public List<double[]> getPolygonCoordinates() {
        return polygonCoordinates;
    }

    public void setPolygonCoordinates(List<double[]> polygonCoordinates) {
        this.polygonCoordinates = polygonCoordinates;
    }

    public int getSatId() {
        return satId;
    }

    public double getReferenceLat() {
        return referenceLat;
    }

    public double getReferenceLon() {
        return referenceLon;
    }

    public double getSurfaceKm2() {
        return surface / 1.0e6;
    }

    public void setSurface(double surface) {
        this.surface = surface;
    }


}
