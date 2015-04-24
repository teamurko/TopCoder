
import java.util.*;

public class SmallPolygons {
    public String[] choosePolygons(int[] coordinates, int maxNumPolygonsAllowed) {
        Point[] points = new Point[coordinates.length / 2];
        for (int i = 0; i < coordinates.length; i += 2) {
            int pointId = i / 2;
            points[pointId] = new Point(pointId, coordinates[i], coordinates[i + 1]);
        }
        List<Polygon> polygons = choosePolygonsImpl(points, maxNumPolygonsAllowed);
        return toJuryOutputFormat(polygons);
    }

    private String[] toJuryOutputFormat(List<Polygon> polygons) {
        String[] res = new String[polygons.size()];
        for (int i = 0; i < res.length; ++i) {
            res[i] = polygons.get(i).serialized();
        }
        return res;
    }

    private List<Polygon> choosePolygonsImpl(Point[] points, int maxPolygons) {
        List<Polygon> res = new ArrayList<>();
        res.add(construct(points, ConstructionStrategy.STAR));
        return res;
    }

    private Polygon construct(Point[] points, ConstructionStrategy strategy) throws IllegalArgumentException {
        if (strategy == ConstructionStrategy.STAR) {
            Polygon convexHull = Polygon.convexHull(points);
            if (convexHull.square() < Utils.epsilon) {
                throw new IllegalArgumentException("Cannot construct 0 area polygon");
            }
            Point[] shiftedPoints = Utils.shiftedBy(points, Utils.center(points));
            Arrays.sort(shiftedPoints, (o1, o2) -> Utils.signum(Math.atan2(o1.y, o1.x) - Math.atan2(o2.y, o2.x)));
            return new Polygon(shiftedPoints);
        } else {
            return null;
        }
    }

    public static void main(String[] args) {
        Scanner in = new Scanner(System.in);
        int numCoordinates = in.nextInt();
        int[] coordinates = new int[numCoordinates];
        for (int i = 0; i < numCoordinates; i++)
            coordinates[i] = in.nextInt();
        SmallPolygons solver = new SmallPolygons();
        String[] serializedPolygons = solver.choosePolygons(coordinates, in.nextInt());
        System.out.println(serializedPolygons.length);
        for (String line : serializedPolygons) {
            System.out.println(line);
        }
        System.out.close();
    }
}

enum ConstructionStrategy {
    STAR
}

class Point {
    final int id;
    final double x;
    final double y;

    Point(int id, double x, double y) {
        this.id = id;
        this.x = x;
        this.y = y;
    }

    Point shifted(double dx, double dy) {
        return new Point(id, x + dx, y + dy);
    }
}

class Polygon {
    final Point[] vertices;

    Polygon(Point...points) {
        this.vertices = points.clone();
    }

    String serialized() {
        return Utils.join(" ", Arrays.stream(vertices).map(p -> p.id).toArray());
    }

    static boolean over(Point a, Point b, Point c) {
        return a.x * (b.y - c.y) + b.x * (c.y - a.y) + c.x * (a.y - b.y) < -Utils.epsilon;
    }

    static boolean under(Point a, Point b, Point c) {
        return a.x * (b.y - c.y) + b.x * (c.y - a.y) + c.x * (a.y - b.y) > Utils.epsilon;
    }

    Point center() {
        double sx = 0;
        double sy = 0;
        for (Point point : vertices) {
            sx += point.x;
            sy += point.y;
        }
        return new Point(Utils.UNDEFINED, sx / vertices.length, sy / vertices.length);
    }

    double square() {
        double sum = 0;
        for (int i = 1; i < vertices.length; i++)
            sum += (vertices[i].x - vertices[i - 1].x) * (vertices[i].y + vertices[i - 1].y);
        sum += (vertices[0].x - vertices[vertices.length - 1].x) * (vertices[0].y + vertices[vertices.length - 1].y);
        return Math.abs(sum) / 2;
    }

    static Polygon convexHull(Point[] points) {
        if (points.length == 1)
            return new Polygon(points);
        Arrays.sort(points, new Comparator<Point>() {
            public int compare(Point o1, Point o2) {
                int value = Double.compare(o1.x, o2.x);
                if (value != 0)
                    return value;
                return Double.compare(o1.y, o2.y);
            }
        });
        Point left = points[0];
        Point right = points[points.length - 1];
        List<Point> up = new ArrayList<>();
        List<Point> down = new ArrayList<>();
        for (Point point : points) {
            if (point == left || point == right || !under(left, point, right)) {
                while (up.size() >= 2 && under(up.get(up.size() - 2), up.get(up.size() - 1), point))
                    up.remove(up.size() - 1);
                up.add(point);
            }
            if (point == left || point == right || !over(left, point, right)) {
                while (down.size() >= 2 && over(down.get(down.size() - 2), down.get(down.size() - 1), point))
                    down.remove(down.size() - 1);
                down.add(point);
            }
        }
        Point[] result = new Point[up.size() + down.size() - 2];
        int index = 0;
        for (Point point : up)
            result[index++] = point;
        for (int i = down.size() - 2; i > 0; i--)
            result[index++] = down.get(i);
        Utils.reverse(result);
        return new Polygon(result);
    }
}

class Utils {
    static final int UNDEFINED = -1;
    static final double epsilon = 1e-7;
    static <T> String join(String delimiter, T[] words) {
       StringBuilder sb = new StringBuilder();
        if (words.length == 0) return "";
        sb.append(words[0]);
        for (int i = 1; i < words.length; ++i) {
            sb.append(delimiter);
            sb.append(words[i]);
        }
        return sb.toString();
    }

    static int signum(double a) {
        if (a == 0.0) return 0;
        return a < 0 ? -1 : 1;
    }

    int vectProductSign(Point o, Point a, Point b) {
        return signum((a.x - o.x) * (b.y - o.y) - (a.y - o.y) * (b.x - o.x));
    }

    static Point[] shiftedBy(Point[] points, Point newOrigin) {
        Point[] res = new Point[points.length];
        for (int i = 0; i < points.length; ++i) {
            res[i] = points[i].shifted(-newOrigin.x, -newOrigin.y);
        }
        return res;
    }

    static void reverse(Point[] points) {
        int l = 0, r = points.length - 1;
        while (l < r) {
            Point t = points[l];
            points[l] = points[r];
            points[r] = t;
            ++l;
            --r;
        }
    }

    static Point center(Point[] points) {
        double sx = 0;
        double sy = 0;
        for (Point point : points) {
            sx += point.x;
            sy += point.y;
        }
        return new Point(Utils.UNDEFINED, sx / points.length, sy / points.length);
    }
}