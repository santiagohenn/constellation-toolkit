package constellation.tools.geometry;

import constellation.tools.math.Pair;

import java.util.List;

public class AAP {

    private long date;
    private int nOfGwsInSight;
    private List<double[]> nonEuclideanCoordinates;
    private List<Double> lats;
    private List<Double> longs;
    List<Integer> gwsInSight;
    private double surfaceInKm2;

    public AAP() {

    }

    public AAP(long date, int nOfGwsInSight, List<Integer> gwsInSight, List<double[]> nonEuclideanCoordinates
            , List<Double> lats, List<Double> longs) {
        this(date, nOfGwsInSight, gwsInSight, nonEuclideanCoordinates, lats, longs, 0);
    }

    public AAP(long date, int nOfGwsInSight, List<Integer> gwsInSight, List<double[]> nonEuclideanCoordinates
            , List<Double> lats, List<Double> longs, double surfaceInKm2) {
        this.date = date;
        this.nOfGwsInSight = nOfGwsInSight;
        this.gwsInSight = gwsInSight;
        this.nonEuclideanCoordinates = nonEuclideanCoordinates;
        this.lats = lats;
        this.longs = longs;
        this.surfaceInKm2 = surfaceInKm2;
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

    public List<double[]> getNonEuclideanCoordinates() {
        return nonEuclideanCoordinates;
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

}
