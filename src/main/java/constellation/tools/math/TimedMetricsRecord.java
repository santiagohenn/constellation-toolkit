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

    public void scale(double factor) {
        for (int i = 0; i < metrics.length; i++) {
            metrics[i] *= factor;
        }
    }

    public long getTime() {
        return time;
    }

    public double[] getMetrics() {
        return metrics;
    }

}
