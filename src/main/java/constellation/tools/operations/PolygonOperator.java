package constellation.tools.operations;

import com.menecats.polybool.Epsilon;
import com.menecats.polybool.PolyBool;
import com.menecats.polybool.models.Polygon;
import satellite.tools.utils.Log;

import java.util.ArrayList;
import java.util.List;

import static com.menecats.polybool.helpers.PolyBoolHelper.epsilon;
import static com.menecats.polybool.helpers.PolyBoolHelper.polygon;

public class PolygonOperator {

    private double epsilonToleranceForOperations;

    public PolygonOperator(double epsilonToleranceForOperations) {
        this.epsilonToleranceForOperations = epsilonToleranceForOperations;
    }

    /**
     * Performs the Intersection of a list of polygons
     *
     * @see <a href="https://github.com/Menecats/polybool-java">Menecats-Polybool</a>
     * @see <a href="https://www.sciencedirect.com/science/article/pii/S0965997813000379">Martinez-Rueda clipping algorithm</a>
     **/
    public Polygon polyIntersect(List<List<double[]>> polygonsToIntersect) {

        if (polygonsToIntersect.stream().allMatch(List::isEmpty)) {
            return new Polygon(new ArrayList<>());
        }

        Polygon intersection = new Polygon();

        int tries = 0;

        while (tries < 3) {

            // Union of all ROI intersections
            try {

                if (polygonsToIntersect.size() < 2) {
                    return intersection;
                }

                Epsilon eps = epsilon(epsilonToleranceForOperations);
                Polygon result = polygon(polygonsToIntersect.get(0));
                PolyBool.Segments segments = PolyBool.segments(eps, result);
                for (int i = 1; i < polygonsToIntersect.size(); i++) {
                    if (polygonsToIntersect.get(i).size() < 3) {
                        continue;
                    }
                    PolyBool.Segments seg2 = PolyBool.segments(eps, polygon(polygonsToIntersect.get(i)));
                    PolyBool.Combined comb = PolyBool.combine(eps, segments, seg2);
                    segments = PolyBool.selectIntersect(comb);
                }

                intersection = PolyBool.polygon(eps, segments);

                if (tries > 0) {
                    Log.warn("Zero-length segment error recovered with epsilon " + epsilonToleranceForOperations);
                }

                break;

            } catch (IndexOutOfBoundsException e1) {
                Log.error("IndexOutOfBoundsException " + e1.getMessage() + " for " + polygonsToIntersect.size() + " regions");
                polygonsToIntersect.forEach(polygon -> Log.error("Poly size " + polygon.size()));
                Log.error(e1.getLocalizedMessage());
            } catch (RuntimeException e2) {
                Log.warn(e2.getMessage());
                Log.warn("RuntimeException. - Increasing epsilon");
                epsilonToleranceForOperations *= 10;
                tries++;
                if (tries == 3) {
                    Log.error("Zero-length segment error could not be recovered.");
                }
            }
        }

        return intersection;

    }


    /**
     * Performs the union of a list of polygons
     *
     * @see <a href="https://github.com/Menecats/polybool-java">Menecats-Polybool</a>
     * @see <a href="https://www.sciencedirect.com/science/article/pii/S0965997813000379">Martinez-Rueda clipping algorithm</a>
     **/
    public Polygon polyUnion(List<List<double[]>> unionQueue) {

        Polygon union = new Polygon();

        if (unionQueue.isEmpty()) {
            return union;
        }  else if (unionQueue.size() == 1) {
            return polygon(unionQueue.get(0));
        }

        int tries = 0;

        while (tries < 3) {

            // Union of all ROI intersections
            try {
                Epsilon eps = epsilon(epsilonToleranceForOperations);
                Polygon result = polygon(unionQueue.get(0));
                PolyBool.Segments segments = PolyBool.segments(eps, result);

                for (int i = 1; i < unionQueue.size(); i++) {
                    PolyBool.Segments seg2 = PolyBool.segments(eps, polygon(unionQueue.get(i)));
                    PolyBool.Combined comb = PolyBool.combine(eps, segments, seg2);
                    segments = PolyBool.selectUnion(comb);
                }

                union = PolyBool.polygon(eps, segments);

                if (tries > 0) {
                    Log.warn("Zero-length segment error recovered with epsilon " + epsilonToleranceForOperations);
                }

                break;

            } catch (IndexOutOfBoundsException e1) {
                Log.error("IndexOutOfBoundsException " + e1.getMessage());
                Log.error("polygons to be or list size: " + unionQueue.size());
            } catch (RuntimeException e2) {
                Log.warn(e2.getMessage());
                Log.warn("RuntimeException. Union size: " + unionQueue.size() + " - Increasing epsilon");
                epsilonToleranceForOperations *= 10;
                tries++;
                if (tries == 3) {
                    Log.error("Zero-length segment error could not be recovered.");
                }
            }
        }

        return union;

    }

    public static boolean pointInPolygon(List<double[]> coordinateList, double[] point) {
        boolean odd = false;

        double[][] polygon = new double[coordinateList.size()][];
        for (int i = 0; i < coordinateList.size(); i++) {
            polygon[i] = coordinateList.get(i);
        }

        //For each edge (In this case for each point of the polygon and the previous one)
        for (int i = 0, j = polygon.length - 1; i < polygon.length; i++) { // Starting with the edge from the last to the first node
            //If a line from the point into infinity crosses this edge
            if (((polygon[i][1] > point[1]) != (polygon[j][1] > point[1])) // One point needs to be above, one below our y coordinate
                    // ...and the edge doesn't cross our Y corrdinate before our x coordinate (but between our x coordinate and infinity)
                    && (point[0] < (polygon[j][0] - polygon[i][0]) * (point[1] - polygon[i][1]) / (polygon[j][1] - polygon[i][1]) + polygon[i][0])) {
                odd = !odd;
            }
            j = i;
        }
        //If the number of crossings was odd, the point is in the polygon
        return odd;
    }

}
