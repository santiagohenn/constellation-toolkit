package constellation.tools.geometry;

import java.util.List;

public class AAP {

    private long date;
    private int nOfGwsInSight;
    private List<double[]> geoCoordinates;
    private List<Double> lats;
    private List<Double> longs;
    List<Integer> gwsInSight;
    private double surfaceInKm2;
    private double referenceLat;
    private double referenceLon;

    public AAP() {

    }

    public AAP(long date, int nOfGwsInSight, List<Integer> gwsInSight, List<double[]> nonEuclideanCoordinates
            , List<Double> lats, List<Double> longs) {
        this(date, nOfGwsInSight, gwsInSight, nonEuclideanCoordinates, lats, longs, 0, 0, 0);
    }

    public AAP(long date, int nOfGwsInSight, List<Integer> gwsInSight, List<double[]> nonEuclideanCoordinates
            , List<Double> lats, List<Double> longs, double referenceLat, double referenceLon) {
        this(date, nOfGwsInSight, gwsInSight, nonEuclideanCoordinates, lats, longs, 0, referenceLat, referenceLon);
    }

    public AAP(long date, int nOfGwsInSight, List<Integer> gwsInSight, List<double[]> nonEuclideanCoordinates
            , List<Double> lats, List<Double> longs, double surfaceInKm2, double referenceLat, double referenceLon) {
        this.date = date;
        this.nOfGwsInSight = nOfGwsInSight;
        this.gwsInSight = gwsInSight;
        this.geoCoordinates = nonEuclideanCoordinates;
        this.lats = lats;
        this.longs = longs;
        this.surfaceInKm2 = surfaceInKm2;
        this.referenceLat = referenceLat;
        this.referenceLon = referenceLon;
    }

    public long getDate() {
        return date;
    }

    public void setDate(long date) {
        this.date = date;
    }

    public int getnOfGwsInSight() {
        return nOfGwsInSight;
    }

    public void setnOfGwsInSight(int nOfGwsInSight) {
        this.nOfGwsInSight = nOfGwsInSight;
    }

    public List<Integer> getGwsInSight() {
        return gwsInSight;
    }

    public void setGwsInSight(List<Integer> gwsInSight) {
        this.gwsInSight = gwsInSight;
    }

    public List<double[]> getGeoCoordinates() {
        return geoCoordinates;
    }

    public List<Double> getLats() {
        return lats;
    }

    public void setLats(List<Double> lats) {
        this.lats = lats;
    }

    public List<Double> getLongs() {
        return longs;
    }

    public void setLongs(List<Double> longs) {
        this.longs = longs;
    }

    public double getSurfaceInKm2() {
        return surfaceInKm2;
    }

    public void setSurfaceInKm2(double surfaceInKm2) {
        this.surfaceInKm2 = surfaceInKm2;
    }

    public void setReference(double referenceLat, double referenceLon) {
        this.referenceLat = referenceLat;
        this.referenceLon = referenceLon;
    }

    public double getReferenceLat() {
        return referenceLat;
    }

    public double getReferenceLon() {
        return referenceLon;
    }

}
