package constellation.tools.geometry;

import constellation.tools.math.Pair;
import constellation.tools.math.Transformations;
import net.sf.geographiclib.*;

import java.awt.geom.Path2D;
import java.util.List;

/**
 * Contains Geographic operations
 * */
public class Geo {

    /**
     * Computes the angular distance over the euclidean plane given two pair of Region objects. Said objects
     * must have reference coordinates, distance will be computed among those coordinates
     *
     * @param r1 the first FOV
     * @param r2 the second FOV
     * @return Double the computed angular distance in degrees
     **/
    public static double computeGeodesic(FOV r1, FOV r2) {
        return computeGeodesic(r1.getReferenceLat(), r1.getReferenceLon(), r2.getReferenceLat(), r2.getReferenceLon());
    }

    /**
     * Computes the angular distance over the euclidean plane given two pair of Coordinates
     *
     * @param lat1 the first FOV's latitude
     * @param lon1 the first FOV's longitude
     * @param lat2 the second FOV's latitude
     * @param lon2 the second FOV's longitude
     * @return a double value for the computed angular distance, in degrees
     **/
    public static double computeGeodesic(double lat1, double lon1, double lat2, double lon2) {
        GeodesicData g = Geodesic.WGS84.Inverse(lat1, lon1, lat2, lon2,
                GeodesicMask.DISTANCE);
        return g.a12;
    }


    /**
     * This method computes the geodetic area of a Path2D.Double object
     *
     * @param polygon The Path2D.Double Object depicting the polygon
     * @return a double value of the computed area in meters squared
     **/
    public static double computeNonEuclideanSurface(Path2D.Double polygon) {
        return computeNonEuclideanSurface(Transformations.polygon2pairList(polygon));
    }

    /**
     * This method computes the geodetic area of a Path2D.Double object
     *
     * @param polygon The Path2D.Double Object depicting the polygon
     * @return a double value of the computed area in kilometers squared
     **/
    @Deprecated
    public static double computeGeodeticAreaKm(Path2D.Double polygon) {
        return computeNonEuclideanSurface(polygon) / 1E6;
    }

    /**
     * This method computes the geodetic area of a list of coordinates, given by Pair objects depicting the polygon's
     * vertices. This method uses the net.sf.geographiclib library.
     *
     * @param pairList A List containing the coordinates of the polygon
     * @return Double the computed area in meters squared
     **/
    public static double computeNonEuclideanSurface(List<Pair> pairList) {

        PolygonArea polygonArea = new PolygonArea(Geodesic.WGS84, false);

        for (Pair pair : pairList) {
            polygonArea.AddPoint(pair.lat, pair.lon);
        }

        PolygonResult result = polygonArea.Compute();
        return Math.abs(result.area);

    }



}
