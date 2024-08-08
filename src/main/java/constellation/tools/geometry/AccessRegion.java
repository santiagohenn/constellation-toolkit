package constellation.tools.geometry;

import java.util.List;

public class AccessRegion {

    private double sspLat;
    private double sspLon;
    private double height;
    private int satId;
    private List<double[]> polygonCoordinates;
    private double surface;
    private double lambdaMax;
    private double x;
    private double y;
    private double z;

    public AccessRegion() {

    }

    public AccessRegion(int satId, double centerLat, double centerLong, List<double[]> polygon) {
        this.satId = satId;
        this.sspLat = centerLat;
        this.sspLon = centerLong;
        this.polygonCoordinates = polygon;
    }

    public AccessRegion(int satId, double centerLat, double centerLong, double height, List<double[]> polygon) {
        this.satId = satId;
        this.sspLat = centerLat;
        this.sspLon = centerLong;
        this.height = height;
        this.polygonCoordinates = polygon;
    }

    public List<double[]> getPolygonCoordinates() {
        return polygonCoordinates;
    }

    public void setPolygonCoordinates(List<double[]> polygonCoordinates) {
        this.polygonCoordinates = polygonCoordinates;
    }

    public double getLambdaMax() {
        return lambdaMax;
    }

    public void setLambdaMax(double lambdaMax) {
        this.lambdaMax = lambdaMax;
    }

    public void setSatPos(double x, double y, double z) {
        this.x = x;
        this.y = y;
        this.z = z;
    }

    public void setSatLLA(double lat, double lon, double height) {
        this.sspLat = lat;
        this.sspLon = lon;
        this.height = height;
    }

    public double getX() {
        return x;
    }

    public double getY() {
        return y;
    }

    public double getZ() {
        return z;
    }

    public int getSatId() {
        return satId;
    }

    public double getSspLat() {
        return sspLat;
    }

    public double getSspLon() {
        return sspLon;
    }

    public double getHeight() {
        return height;
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
