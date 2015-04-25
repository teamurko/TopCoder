
import java.util.*;
import java.util.stream.Collectors;

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
        List<Point[]> groups = split(points, maxPolygons);
        res.addAll(groups.stream().map(group -> construct(group, ConstructionStrategy.STAR)).collect(Collectors.toList()));
        return res;
    }

    private List<Point[]> split(Point[] points, int maxGroups) {
        double[][] dist = new double[points.length][points.length];
        for (int i = 0; i < points.length; ++i) {
            for (int j = i + 1; j < points.length; ++j) {
                int iId = points[i].id;
                int jId = points[j].id;
                dist[iId][jId] = dist[jId][iId] = points[i].distanceTo(points[j]);
            }
        }
        Polygon convexHull = Polygon.convexHull(points);
        //TODO use better initial centers
        maxGroups = Math.min(maxGroups, convexHull.vertices.length);
        Point[] centers = choose(convexHull.vertices.clone(), maxGroups);
        boolean updated = true;
        for (int iter = 0; iter < 100 && updated; ++iter) {
            List<Point>[] groups = distribute(points, centers);
            updated = false;
            for (int i = 0; i < maxGroups; ++i) {
                Point center = Utils.center(groups[i]);
                if (centers[i].distanceTo(center) > Utils.epsilon) {
                    updated = true;
                    centers[i] = center;
                }
            }
        }

        List<Point> newCenters = new ArrayList<>();
        for (List<Point> group : distribute(points, centers)) {
            if (group.size() > 3) {
                newCenters.add(Utils.center(group));
            }
        }
        List<Point[]> res = new ArrayList<>();
        centers = new Point[newCenters.size()];
        newCenters.toArray(centers);
        for (List<Point> group : distribute(points, centers)) {
            Point[] pts = new Point[group.size()];
            group.toArray(pts);
            res.add(pts);
        }
        return res;
    }

    private static List<Point>[] distribute(Point[] points, Point[] centers) {
        List<Point>[] res = new ArrayList[centers.length];
        for (int i = 0; i < res.length; ++i) {
            res[i] = new ArrayList<>();
        }
        for (Point p : points) {
            res[Utils.closest(centers, p)].add(p);
        }
        return res;
    }

    private static Point[] choose(Point[] points, int numPoints) {
        Point[] res = new Point[numPoints];
        int ptr = 0;
        while (numPoints > 0) {
            int i = Utils.random.nextInt(numPoints);
            res[ptr++] = points[i];
            Point t = points[numPoints - 1];
            points[numPoints - 1] = points[i];
            points[i] = t;
            --numPoints;
        }
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

    static TestData readData(boolean generate, Scanner in) {
        if (generate) {
            return Utils.generate(in.next());
        } else {
            int numCoordinates = in.nextInt();
            int[] coordinates = new int[numCoordinates];
            for (int i = 0; i < numCoordinates; i++) {
                coordinates[i] = in.nextInt();
            }
            int maxPolygons = in.nextInt();
            return new TestData(coordinates, maxPolygons);
        }
    }

    public static void main(String[] args) {
        Scanner in = new Scanner(System.in);
        boolean generate = args.length > 0 && "generate".equals(args[0]);
        TestData testData = readData(generate, in);
        SmallPolygons solver = new SmallPolygons();
        String[] serializedPolygons = solver.choosePolygons(testData.coordinates, testData.maxPolygons);
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

class TestData {
    final int[] coordinates;
    final int maxPolygons;

    TestData(int[] coordinates, int maxPolygons) {
        this.coordinates = coordinates;
        this.maxPolygons = maxPolygons;
    }
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

    double distanceTo(Point p) {
        return Math.hypot(p.x - x, p.y - y);
    }
}

class DisjointSetUnion {
    final int[] parent;
    final Random random = new Random();
    DisjointSetUnion(int size) {
        parent = new int[size];
        for (int i = 0; i < size; ++i) parent[i] = i;
    }

    int getParent(int v) {
        if (parent[v] != v) parent[v] = getParent(parent[v]);
        return parent[v];
    }

    void merge(int u, int v) {
        u = getParent(u);
        v = getParent(v);
        if (u != v) {
            if (random.nextBoolean()) {
                parent[u] = v;
            } else {
                parent[v] = u;
            }
        }
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
        return new Point(Utils.undefined, sx / vertices.length, sy / vertices.length);
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
    static final int undefined = -1;
    static final double epsilon = 1e-7;
    static final Random random = new Random();

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

    static int closest(Point[] points, Point target) {
        int res = 0;
        double dst = target.distanceTo(points[0]);
        for (int i = 1; i < points.length; ++i) {
            double candidate = target.distanceTo(points[i]);
            if (candidate < dst) {
                dst = candidate;
                res = i;
            }
        }
        return res;
    }

    static int vectProductSign(Point o, Point a, Point b) {
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
        return new Point(Utils.undefined, sx / points.length, sy / points.length);
    }

    static Point center(List<Point> points) {
        double sx = 0;
        double sy = 0;
        for (Point point : points) {
            sx += point.x;
            sy += point.y;
        }
        return new Point(Utils.undefined, sx / points.size(), sy / points.size());
    }

    static Point[] coord2Points(int[] coords) {
        Point[] points = new Point[coords.length / 2];
        for (int i = 0; i < coords.length; i += 2) {
            points[i] = new Point(i / 2, coords[i], coords[i + 1]);
        }
        return points;
    }

    static TestData generate(String seed) {
        SmallPolygonsVis spv = new SmallPolygonsVis();
        spv.generate(seed);
        return new TestData(spv.pointsPar, spv.N);
    }
}