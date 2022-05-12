package constellation.tools.math;

public class Pair {

    public double lat;
    public double lon;

    public Pair(double first, double second) {
        this.lat = first;
        this.lon = second;
    }

    public Pair() {

    }

    public double[] getPoint() {
        return new double[]{lat, lon};
    }

    @Override
    public String toString() {
        return lat + "," + lon;
    }
}