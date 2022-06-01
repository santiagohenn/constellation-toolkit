package constellation.tools.geometry;

import java.awt.geom.Path2D;
import java.util.List;

public class FOV {

    private final double referenceLat;
    private final double referenceLon;
    private final int satId;
    private final Path2D.Double polygon;
    private List<double[]> polygonCoordinates;
    private double surface;

    public FOV(int satId, double centerLat, double centerLong, Path2D.Double polygon) {
        this.satId = satId;
        this.referenceLat = centerLat;
        this.referenceLon = centerLong;
        this.polygon = polygon;
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

    public Path2D.Double getPolygon() {
        return polygon;
    }

    public double getSurface() {
        return surface;
    }

    public double getSurfaceKm2() {
        return surface / 1.0e6;
    }

    public void setSurface(double surface) {
        this.surface = surface;
    }


}
