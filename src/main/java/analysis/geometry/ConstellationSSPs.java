package analysis.geometry;

import analysis.math.Pair;

import java.util.ArrayList;
import java.util.List;

public class ConstellationSSPs {

    private final long date;
    private final List<Double> lats = new ArrayList<>();
    private final List<Double> longs = new ArrayList<>();
    private List<Pair> sspCoordinates = new ArrayList<>();

    public ConstellationSSPs(long date) {
        this.date = date;
    }

    public ConstellationSSPs(long date, List<Pair> sspCoordinates) {
        this.date = date;
        this.sspCoordinates = sspCoordinates;
        sspCoordinates.forEach(pair -> {
            lats.add(pair.lat);
            longs.add(pair.lon);
        });
    }

    public void addSSP(Pair pair) {
        this.sspCoordinates.add(pair);
        lats.add(pair.lat);
        longs.add(pair.lon);
    }

    public void addSSP(int id, Pair coordinates) {
        this.sspCoordinates.add(id, coordinates);
    }

    public long getDate() {
        return date;
    }

    public List<Double> getLats() {
        return lats;
    }

    public List<Double> getLongs() {
        return longs;
    }
}
