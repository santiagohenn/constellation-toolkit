package analysis.geometry;

import java.awt.geom.Path2D;
import java.awt.geom.PathIterator;
import java.util.ArrayList;
import java.util.List;

public class FOV {

    private final double referenceLat;
    private final double referenceLon;
    private final int satId;
    private final Path2D.Double polygon;
    private double surface;

    public FOV(int satId, double centerLat, double centerLong, Path2D.Double polygon) {
        this.satId = satId;
        this.referenceLat = centerLat;
        this.referenceLon = centerLong;
        this.polygon = polygon;
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

    // TODO: Remove if unused
    public List<Pair> getPolygonAsPairList() {

        List<Pair> pairList = new ArrayList<>();

        PathIterator iterator = polygon.getPathIterator(null);

        final double[] point = new double[2];

        while (!iterator.isDone()) {

            int type = iterator.currentSegment(point);

            if (type == PathIterator.SEG_CLOSE) {
                break;
            }

            pairList.add(new Pair(point[0], point[1]));
            iterator.next();

        }

        return pairList;

    }

}
