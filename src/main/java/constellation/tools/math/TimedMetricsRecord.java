package constellation.tools.math;

public class TimedMetricsRecord {

    long time;
    double[] metrics;

    public TimedMetricsRecord(long time, int metricsCount) {
        this.time = time;
        metrics = new double[metricsCount];
    }

    public void addMetric(int assetID, double surfaceValue) {
        metrics[assetID] = surfaceValue;
    }

    public long getTime() {
        return time;
    }

    public double[] getMetrics() {
        return metrics;
    }

}
