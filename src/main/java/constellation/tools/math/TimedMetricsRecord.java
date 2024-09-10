package constellation.tools.math;

import java.util.Arrays;

public class TimedMetricsRecord {

    long time;
    double[] metrics;

    public TimedMetricsRecord(long time, int metricsCount) {
        this.time = time;
        metrics = new double[metricsCount];
        Arrays.fill(metrics, 0.0);
    }

    public void addMetric(int index, double value) {
        metrics[index] = value;
    }

    public void accumulateMetric(int index, double value) {
        metrics[index] += value;
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

    @Override
    public String toString() {
        return time + ", " + Arrays.toString(metrics);
    }

}
