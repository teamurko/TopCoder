import javafx.util.Pair;

import javax.swing.*;
import java.awt.*;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.WindowEvent;
import java.awt.event.WindowListener;
import java.awt.image.BufferedImage;
import java.io.*;
import java.security.SecureRandom;
import java.util.*;
import java.util.List;
import java.util.stream.Collectors;

public class SmallPolygons {
    private final static SimplePolygonizationBuilder polyBuilder = new SimplePolygonizationBuilder();
    private final Set<ConstructionStrategy> strategies;

    public SmallPolygons() {
        this(EnumSet.allOf(ConstructionStrategy.class));
    }

    public SmallPolygons(Set<ConstructionStrategy> strategies) {
        this.strategies = strategies;
    }

    public String[] choosePolygons(int[] coordinates, int maxNumPolygonsAllowed) {
        Utils.timer = new Timer();
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
        res.addAll(groups.stream().map(group -> polyBuilder.build(group, strategies)).collect(Collectors.toList()));
        return res;
    }

    private List<Point[]> split(Point[] points, int maxGroups) {
        maxGroups = Math.min(maxGroups, points.length / 3);
        Point[] centers = Utils.choose(points, maxGroups);
        boolean updated = true;
        for (int iter = 0; iter < 200 && updated; ++iter) {
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
            if (group.size() > 2) {
                newCenters.add(Utils.center(group));
            }
        }
        List<Point[]> res = new ArrayList<>();
        centers = new Point[newCenters.size()];
        newCenters.toArray(centers);
        List<Point>[] groups = distribute(points, centers);
        maxGroups -= groups.length;
        PriorityQueue<List<Point>> groupsBySize = new PriorityQueue<>((o1, o2) -> -o1.size() + o2.size());
        for (List<Point> group : groups) {
            groupsBySize.add(group);
        }
        while (maxGroups > 0 && groupsBySize.peek().size() >= 6) {
            List<Point> group = groupsBySize.poll();
            List<Point> first = new ArrayList<>();
            List<Point> second = new ArrayList<>();
            if (trySplit(group, first, second)) {
                groupsBySize.add(first);
                groupsBySize.add(second);
            }
            --maxGroups;
        }
        for (List<Point> group : groupsBySize) {
            Point[] pts = new Point[group.size()];
            group.toArray(pts);
            res.add(pts);
        }
        return res;
    }

    static boolean trySplit(List<Point> group, List<Point> first, List<Point> second) {
        Collections.sort(group, (o1, o2) -> {
            int res = Utils.signum(o1.x - o2.x);
            if (res != 0) return res;
            return Utils.signum(o1.y - o2.y);
        });
        for (Point p : group) {
            if (first.size() * 2 < group.size()) {
                first.add(p);
            } else {
                second.add(p);
            }
        }
        return first.size() >= 3 && second.size() >= 3;
    }

    static List<Point>[] distribute(Point[] points, Point[] centers) {
        List<Point>[] res = new ArrayList[centers.length];
        for (int i = 0; i < res.length; ++i) {
            res[i] = new ArrayList<>();
        }
        for (Point p : points) {
            res[Utils.closest(centers, p)].add(p);
        }
        return res;
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
        Set<ConstructionStrategy> strategies;
        if (args.length > 0 && "improve".equals(args[0])) {
            strategies = EnumSet.allOf(ConstructionStrategy.class);
        } else {
            strategies = EnumSet.of(
                    ConstructionStrategy.STAR,
                    ConstructionStrategy.STAR_SKEWED,
                    ConstructionStrategy.TRY_BRUTE_FORCE,
                    ConstructionStrategy.A1,
                    ConstructionStrategy.A2,
                    ConstructionStrategy.A3,
                    ConstructionStrategy.A4);
        }
        boolean generate = args.length > 1 && "generate".equals(args[1]);
        Scanner in = new Scanner(System.in);
        TestData testData = readData(generate, in);
        SmallPolygons solver = new SmallPolygons(strategies);
        String[] serializedPolygons = solver.choosePolygons(testData.coordinates, testData.maxPolygons);
        System.out.println(serializedPolygons.length);
        for (String line : serializedPolygons) {
            System.out.println(line);
        }
        System.out.close();
    }
}

class BruteForcer {
    Polygon run(Point[] points) {
        PolygonHolder polygonHolder = new PolygonHolder();
        Point[] res = new Point[points.length];
        res[0] = points[0];
        boolean[] taken = new boolean[points.length];
        taken[0] = true;
        runImpl(res, 1, points, taken, polygonHolder);
        return polygonHolder.getPolygon();
    }

    private void runImpl(Point[] res, int i, Point[] points, boolean[] taken, PolygonHolder polygonHolder) {
        if (i == res.length) {
            boolean ok = true;
            for (int k = 1; k < i; ++k) {
                if (Utils.intersect(res[k - 1], res[k], res[i - 1], res[0])) {
                    ok = false;
                    break;
                }
            }
            if (ok) {
                polygonHolder.update(new Polygon(res));
            }
            return;
        }
        for (int j = 0; j < points.length; ++j) {
            if (!taken[j]) {
                boolean ok = true;
                for (int k = 1; k < i; ++k) {
                    if (Utils.intersect(res[k - 1], res[k], res[i - 1], points[j])) {
                        ok = false;
                        break;
                    }
                }
                if (ok) {
                    taken[j] = true;
                    res[i] = points[j];
                    runImpl(res, i + 1, points, taken, polygonHolder);
                    taken[j] = false;
                }
            }
        }
    }
}

class PolygonHolder {
    private Polygon best;
    private double bestArea = Utils.infinity;
    void update(Polygon candidate) {
        if (candidate != null) {
            double area = candidate.area();
            if (area < bestArea) {
                bestArea = area;
                best = candidate;
            }
        }
    }
    double getArea() { return bestArea; }
    Polygon getPolygon() { return best; }
}

interface PolygonConstructor {
    Polygon build(Point[] points);
}

class SimplePolygonizationBuilder {
    private static final Map<ConstructionStrategy, PolygonConstructor> strategyToConstructor = new HashMap<>();
    static {
        Utils.addConstructors(strategyToConstructor);
    }
    Polygon build(Point[] points, Set<ConstructionStrategy> strategies) throws IllegalArgumentException {
        Polygon convexHull = Polygon.convexHull(points);
        if (convexHull.area() < Utils.epsilon) {
            throw new IllegalArgumentException("Cannot construct 0 area polygon");
        }
        PolygonHolder polygonHolder = new PolygonHolder();
        strategyToConstructor.entrySet().stream().filter(
                entry -> strategies.contains(entry.getKey())).forEach(
                    entry -> polygonHolder.update(entry.getValue().build(points)));
        return polygonHolder.getPolygon();
    }
}

enum ConstructionStrategy {
    STAR,
    STAR_SKEWED,
    TRY_BRUTE_FORCE,
    A1,
    A2,
    A3,
    A4,
    A5,
    A6
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

    double distanceTo(Point a, Point b) {
        double A = -(a.x - x) * (b.x - a.x) - (a.y - y) * (b.y - a.y);
        double t = A / ((b.x - a.x) * (b.x - a.x) + (b.y - a.y) * (b.y - a.y));
        if (0 < t && t < 1) {
            return distanceTo(new Point(Utils.undefined, a.x + (b.x - a.x) * t, a.y + (b.y - a.y) * t));
        }
        return Math.min(a.distanceTo(this), b.distanceTo(this));
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

    Polygon(List<Point> points) {
        vertices = new Point[points.size()];
        points.toArray(vertices);
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

    double area() {
        double sum = 0;
        for (int i = 1; i < vertices.length; i++)
            sum += (vertices[i].x - vertices[i - 1].x) * (vertices[i].y + vertices[i - 1].y);
        sum += (vertices[0].x - vertices[vertices.length - 1].x) * (vertices[0].y + vertices[vertices.length - 1].y);
        return Math.abs(sum) / 2;
    }

    boolean contains(Point point, boolean strict) {
        for (int i = 0; i < vertices.length; ++i) {
            int j = i + 1;
            if (j == vertices.length) j = 0;
            int sign = Utils.vectProductSign(point, vertices[i], vertices[j]);
            if (sign < 0) return false;
            if (sign == 0 && strict) return false;
        }
        return true;
    }

    static Polygon convexHull(List<Point> points) {
        Point[] pts = new Point[points.size()];
        points.toArray(pts);
        return convexHull(pts);
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

    Polygon extendConvexHull(Point p) {
        List<Point> res = new ArrayList<>();
        for (int i = 0; i < vertices.length; ++i) {
            int j = i - 1;
            if (j < 0) j = vertices.length - 1;
            int k = i + 1;
            if (k == vertices.length) k = 0;
            if (Utils.vectProductSign(p, vertices[j], vertices[i]) >= 0 || Utils.vectProductSign(p, vertices[i], vertices[k]) >= 0) {
                res.add(vertices[i]);
            }
        }
        for (int i = 0; i < res.size(); ++i) {
            int j = i + 1;
            if (j == res.size()) j = 0;
            if (Utils.vectProductSign(p, vertices[i], vertices[j]) < 0) {
                res.add(j, p);
                break;
            }
        }
        return new Polygon(res);
    }

    boolean containsAny(Point[] points) {
        for (Point p : points) {
            if (contains(p, true)) {
                return true;
            }
        }
        return false;
    }
}

class Timer {
    private static final long timeLimit = 10000; //ms
    private static final long waitTimeLimit = 2000; //ms
    private long startTime;
    Timer() { startTime = System.nanoTime(); }
    boolean hasTime() {
        return (System.nanoTime() - startTime) / 1000000 < timeLimit - waitTimeLimit;
    }
}
class Utils {
    static final int undefined = -1;
    static final double epsilon = 1e-7;
    static final Random random = new Random(123);
    static final double infinity = 1e10;
    static final int brute_force_limit = 10;
    static int random_search_limit = 500;

    private static final BruteForcer bruteForcer = new BruteForcer();
    static Timer timer;

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

    static boolean intersect(Point a, Point b, Point c, Point d) {
        //(b.x - a.x) * t1 + (c.x - d.x) * t2 = c.x - a.x
        if (Math.max(a.x, b.x) + epsilon < Math.min(c.x, d.x)) return false;
        if (Math.max(c.x, d.x) + epsilon < Math.min(a.x, b.x)) return false;
        if (Math.max(a.y, b.y) + epsilon < Math.min(c.y, d.y)) return false;
        if (Math.max(c.y, d.y) + epsilon < Math.min(a.y, b.y)) return false;
        double a11 = b.x - a.x;
        double a12 = c.x - d.x;
        double b1 = c.x - a.x;
        double a21 = b.y - a.y;
        double a22 = c.y - d.y;
        double b2 = c.y - a.y;
        double determinant = det(a11, a12, a21, a22);
        if (Math.abs(determinant) < epsilon) {
            return vectProductSign(a, b, c) == 0;
        }
        double t1 = det(b1, a12, b2, a22) / determinant;
        double t2 = det(a11, b1, a21, b2) / determinant;
        if (t1 < -epsilon || t1 > 1 + epsilon || t2 < -epsilon || t2 > 1 + epsilon) return false;
        return (epsilon < t1 && t1 < 1 - epsilon) || (epsilon < t2 && t2 < 1 - epsilon);
    }

    static double det(double a11, double a12, double a21, double a22) {
        return a11 * a22 - a12 * a21;
    }

    static double canonicalAngle(double angle) {
        while (angle > Math.PI)
            angle -= 2 * Math.PI;
        while (angle < -Math.PI)
            angle += 2 * Math.PI;
        return angle;
    }

    static Point[] choose(Point[] pts, int numPoints) {
        Point[] points = pts.clone();
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

    static Point[] randomTriangle(Point[] points) {
        if (points.length < 3) {
            throw new IllegalArgumentException("Not enough points to construct triangle");
        }
        Point[] res = choose(points, 3);
        if (vectProductSign(res[0], res[1], res[2]) < 0) {
            reverse(res);
        }
        return res;
    }

    static Point[] randomGreedyTriangle(Point[] points) {
        if (points.length < 3) {
            throw new IllegalArgumentException("Not enough points to construct triangle");
        }
        Point o = points[Utils.random.nextInt(points.length)];
        Point a = null;
        Point b = null;
        for (Point p : points) {
            if (p != o) {
                if (a == null) {
                    a = p;
                    continue;
                }
                if (a.distanceTo(o) > p.distanceTo(o)) {
                    b = a;
                    a = p;
                } else {
                    if (b == null || b.distanceTo(o) > p.distanceTo(o)) {
                        b = p;
                    }
                }
            }
        }
        if (vectProductSign(o, a, b) < 0) {
            return new Point[]{o, b, a};
        }
        return new Point[]{o, a, b};
    }

    static Point randomPoint(Point[] points) {
        return points[Utils.random.nextInt(points.length)];
    }

    static double area(Point a, Point b, Point c) {
        return Math.abs(det(b.x - a.x, b.y - a.y, c.x - a.x, c.y - a.y));
    }

    static Point[] feasiblePoints(List<Point> vertices, Point[] points) {
        Polygon convexHull = Polygon.convexHull(vertices);
        List<Point> res = new ArrayList<>();
        for (Point point : points) {
            Polygon ech = convexHull.extendConvexHull(point);
            boolean contains = false;
            for (Point p : points) {
                if (p == point) continue;
                if (ech.contains(p, true)) {
                    contains = true;
                }
            }
            if (!contains) {
                res.add(point);
            }
        }
        Point[] pts = new Point[res.size()];
        res.toArray(pts);
        return pts;
    }

    static List<Integer> visibleEdges(List<Point> vertices, Point point) {
        List<Integer> res = new ArrayList<>();
        boolean[] visible = new boolean[vertices.size()];
        for (int i = 0; i < vertices.size(); ++i) {
            if (!intersectAny(point, vertices.get(i), vertices)) {
                visible[i] = true;
            }
        }
        for (int i = 0; i < vertices.size(); ++i) {
            int j = i + 1;
            if (j == vertices.size()) j = 0;
            if (visible[i] && visible[j]) {
                res.add(i);
            }
        }
        return res;
    }

    static boolean intersectAny(Point a, Point b, List<Point> vertices) {
        if (intersect(a, b, vertices.get(0), vertices.get(vertices.size() - 1))) return true;
        for (int i = 0; i + 1 < vertices.size(); ++i) {
            if (intersect(a, b, vertices.get(i), vertices.get(i + 1))) return true;
        }
        return false;
    }

    static int greedyEdge(List<Point> vertices, Point point) {
        int res = Utils.undefined;
        double area = Utils.infinity;
        for (int i : visibleEdges(vertices, point)) {
            int j = i + 1;
            if (j == vertices.size()) j = 0;
            double a = area(point, vertices.get(i), vertices.get(j));
            if (a > epsilon && area > a) {
                area = a;
                res = i;
            }
        }
        return res;
    }

    static Pair<Integer, Integer> greedyArea(List<Point> vertices, Point[] feasiblePoints) {
        int resSegment = Utils.undefined;
        int resPoint = Utils.undefined;
        double area = Utils.infinity;
        for (int i = 0; i < feasiblePoints.length; ++i) {
            int candidate = greedyEdge(vertices, feasiblePoints[i]);
            if (candidate == Utils.undefined) continue;
            double a = area(feasiblePoints[i], vertices.get(candidate), vertices.get((candidate + 1) % vertices.size()));
            if (a > epsilon && a < area) {
                area = a;
                resPoint = i;
                resSegment = candidate;
            }
        }
        return new Pair<>(resSegment, resPoint);
    }

    static double distanceToPolygon(List<Point> vertices, Point point) {
        double res = Utils.infinity;
        for (int i = 0; i < vertices.size(); ++i) {
            res = Math.min(res, point.distanceTo(vertices.get(i), vertices.get((i + 1) % vertices.size())));
        }
        return res;
    }

    static Point closest(Point[] points, List<Point> vertices) {
        double dist = Utils.infinity;
        Point res = null;
        for (Point p : points) {
            double dst = distanceToPolygon(vertices, p);
            if (res == null || dst < dist) {
                dist = dst;
                res = p;
            }
        }
        return res;
    }

    static void addConstructors(Map<ConstructionStrategy, PolygonConstructor> strategyToConstructor) {
        strategyToConstructor.put(ConstructionStrategy.STAR, points -> {
            if (!Utils.timer.hasTime()) return null;
            Point[] shiftedPoints = Utils.shiftedBy(points, Utils.center(points));
            Arrays.sort(shiftedPoints, (o1, o2) -> Utils.signum(Math.atan2(o1.y, o1.x) - Math.atan2(o2.y, o2.x)));
            return new Polygon(shiftedPoints);
        });
        strategyToConstructor.put(ConstructionStrategy.STAR_SKEWED, points -> {
            if (!Utils.timer.hasTime()) return null;
            PolygonHolder polygonHolder = new PolygonHolder();
            Point center = Utils.center(points);
            for (Point candidate : points) {
                double length = Math.hypot(center.x - candidate.x, center.y - candidate.y);
                Point skewedCenter = candidate.shifted((center.x - candidate.x) / length, (center.y - candidate.y) / length);
                Point[] shiftedPoints = Utils.shiftedBy(points, skewedCenter);
                Arrays.sort(shiftedPoints, (o1, o2) -> Utils.signum(Math.atan2(o1.y, o1.x) - Math.atan2(o2.y, o2.x)));
                polygonHolder.update(new Polygon(shiftedPoints));
            }
            return polygonHolder.getPolygon();
        });
        strategyToConstructor.put(ConstructionStrategy.TRY_BRUTE_FORCE, points -> {
            if (points.length > Utils.brute_force_limit || !Utils.timer.hasTime()) return null;
            return bruteForcer.run(points);
        });
        strategyToConstructor.put(ConstructionStrategy.A1, points -> {
            boolean pass = points.length > Utils.random_search_limit || !Utils.timer.hasTime();
            if (pass) return null;
            Point[] initialTriangle = randomGreedyTriangle(points);
            List<Point> vertices = new ArrayList<>(3);
            for (Point p : initialTriangle) vertices.add(p);
            Point[] pts = new Point[points.length - 3];
            int i = 0;
            for (Point p : points) {
                if (!vertices.contains(p)) {
                    pts[i++] = p;
                }
            }

            while (pts.length > 0) {
                Point[] feasiblePoints = feasiblePoints(vertices, pts);
                Point p = closest(feasiblePoints, vertices);
                if (p == null) return null;
                int edge = greedyEdge(vertices, p) + 1;
                if (edge == vertices.size()) edge = 0;
                vertices.add(edge, p);
                for (i = 0; i + 1 < pts.length; ++i) {
                    if (pts[i] == p) {
                        pts[i] = pts[pts.length - 1];
                        break;
                    }
                }
                pts = Arrays.copyOf(pts, pts.length - 1);
            }
            if (!valid(vertices)) return null;
            return new Polygon(vertices);
        });
        strategyToConstructor.put(ConstructionStrategy.A2, points -> {
            boolean pass = points.length > Utils.random_search_limit || !Utils.timer.hasTime();
            if (pass) return null;
            Point[] initialTriangle = randomGreedyTriangle(points);
            List<Point> vertices = new ArrayList<>(3);
            for (Point p : initialTriangle) vertices.add(p);
            Point[] pts = new Point[points.length - 3];
            int i = 0;
            for (Point p : points) {
                if (!vertices.contains(p)) {
                    pts[i++] = p;
                }
            }

            while (pts.length > 0) {
                Point[] feasiblePoints = feasiblePoints(vertices, pts);
                if (feasiblePoints.length == 0) return null;
                Point p = randomPoint(feasiblePoints);
                if (p == null) return null;
                int edge = greedyEdge(vertices, p) + 1;
                if (edge == vertices.size()) edge = 0;
                vertices.add(edge, p);
                for (i = 0; i + 1 < pts.length; ++i) {
                    if (pts[i] == p) {
                        pts[i] = pts[pts.length - 1];
                        break;
                    }
                }
                pts = Arrays.copyOf(pts, pts.length - 1);
            }
            if (!valid(vertices)) return null;
            return new Polygon(vertices);
        });
        strategyToConstructor.put(ConstructionStrategy.A3, points -> {
            boolean pass = points.length > Utils.random_search_limit || !Utils.timer.hasTime();
            if (pass) return null;
            Point[] initialTriangle = randomGreedyTriangle(points);
            List<Point> vertices = new ArrayList<>(3);
            for (Point p : initialTriangle) vertices.add(p);
            Point[] pts = new Point[points.length - 3];
            int i = 0;
            for (Point p : points) {
                if (!vertices.contains(p)) {
                    pts[i++] = p;
                }
            }

            while (pts.length > 0) {
                Point[] feasiblePoints = feasiblePoints(vertices, pts);
                if (feasiblePoints.length == 0) return null;
                Pair<Integer, Integer> t = greedyArea(vertices, feasiblePoints);
                if (t.getValue() == Utils.undefined || t.getKey() == Utils.undefined) return null;
                Point p = feasiblePoints[t.getValue()];
                int edge = t.getKey() + 1;
                if (edge == vertices.size()) edge = 0;
                vertices.add(edge, p);
                for (i = 0; i + 1 < pts.length; ++i) {
                    if (pts[i] == p) {
                        pts[i] = pts[pts.length - 1];
                        break;
                    }
                }
                pts = Arrays.copyOf(pts, pts.length - 1);
            }
            if (!valid(vertices)) return null;
            return new Polygon(vertices);
        });
        strategyToConstructor.put(ConstructionStrategy.A4, points -> {
            boolean pass = points.length > Utils.random_search_limit || !Utils.timer.hasTime();
            if (pass) return null;
            Point[] initialTriangle = randomTriangle(points);
            List<Point> vertices = new ArrayList<>(3);
            for (Point p : initialTriangle) vertices.add(p);
            Point[] pts = new Point[points.length - 3];
            int i = 0;
            for (Point p : points) {
                if (!vertices.contains(p)) {
                    pts[i++] = p;
                }
            }

            while (pts.length > 0) {
                Point[] feasiblePoints = feasiblePoints(vertices, pts);
                if (feasiblePoints.length == 0) return null;
                Point p = closest(feasiblePoints, vertices);
                if (p == null) return null;
                int edge = greedyEdge(vertices, p) + 1;
                if (edge == vertices.size()) edge = 0;
                vertices.add(edge, p);
                for (i = 0; i + 1 < pts.length; ++i) {
                    if (pts[i] == p) {
                        pts[i] = pts[pts.length - 1];
                        break;
                    }
                }
                pts = Arrays.copyOf(pts, pts.length - 1);
            }
            if (!valid(vertices)) return null;
            return new Polygon(vertices);
        });
        strategyToConstructor.put(ConstructionStrategy.A5, points -> {
            boolean pass = points.length > Utils.random_search_limit || !Utils.timer.hasTime();
            if (pass) return null;
            Point[] initialTriangle = randomTriangle(points);
            List<Point> vertices = new ArrayList<>(3);
            for (Point p : initialTriangle) vertices.add(p);
            Point[] pts = new Point[points.length - 3];
            int i = 0;
            for (Point p : points) {
                if (!vertices.contains(p)) {
                    pts[i++] = p;
                }
            }

            while (pts.length > 0) {
                Point[] feasiblePoints = feasiblePoints(vertices, pts);
                if (feasiblePoints.length == 0) return null;
                Point p = randomPoint(feasiblePoints);
                if (p == null) return null;
                int edge = greedyEdge(vertices, p) + 1;
                if (edge == vertices.size()) edge = 0;
                vertices.add(edge, p);
                for (i = 0; i + 1 < pts.length; ++i) {
                    if (pts[i] == p) {
                        pts[i] = pts[pts.length - 1];
                        break;
                    }
                }
                pts = Arrays.copyOf(pts, pts.length - 1);
            }
            if (!valid(vertices)) return null;
            return new Polygon(vertices);
        });
        strategyToConstructor.put(ConstructionStrategy.A6, points -> {
            boolean pass = points.length > Utils.random_search_limit || !Utils.timer.hasTime();
            if (pass) return null;
            Point[] initialTriangle = randomTriangle(points);
            List<Point> vertices = new ArrayList<>(3);
            for (Point p : initialTriangle) vertices.add(p);
            Point[] pts = new Point[points.length - 3];
            int i = 0;
            for (Point p : points) {
                if (!vertices.contains(p)) {
                    pts[i++] = p;
                }
            }

            while (pts.length > 0) {
                Point[] feasiblePoints = feasiblePoints(vertices, pts);
                if (feasiblePoints.length == 0) return null;
                Pair<Integer, Integer> t = greedyArea(vertices, feasiblePoints);
                if (t.getValue() == Utils.undefined || t.getKey() == Utils.undefined) return null;
                Point p = feasiblePoints[t.getValue()];
                int edge = t.getKey() + 1;
                if (edge == vertices.size()) edge = 0;
                vertices.add(edge, p);
                for (i = 0; i + 1 < pts.length; ++i) {
                    if (pts[i] == p) {
                        pts[i] = pts[pts.length - 1];
                        break;
                    }
                }
                pts = Arrays.copyOf(pts, pts.length - 1);
            }
            if (!valid(vertices)) return null;
            return new Polygon(vertices);
        });
    }

    static boolean valid(List<Point> vertices) {
        for (int i = 0; i + 1 < vertices.size(); ++i) {
            for (int j = i + 1; j < vertices.size(); ++j) {
                int k = j + 1;
                if (k == vertices.size()) k = 0;
                if (intersect(vertices.get(i), vertices.get(i + 1), vertices.get(j), vertices.get(k))) return false;
            }
        }
        return true;
    }

    static boolean threePointsOnLine(Point[] points) {
        for (int i = 0; i < points.length; ++i) {
            for (int j = i + 1; j < points.length; ++j) {
                for (int k = j + 1; k < points.length; ++k) {
                    if (vectProductSign(points[i], points[j], points[k]) == 0) return true;
                }
            }
        }
        return false;
    }
}

// ------------- class Point ------------------------------
class Pnt {
    public int x,y;
    public Pnt() {};
    public Pnt(int x1, int y1) {
        x = x1;
        y = y1;
    }
    public boolean equals(Pnt other) {
        return (x == other.x && y == other.y);
    }
}

// ------------- class G2D --------------------------------
class G2D {
    public static Pnt substr(Pnt p1, Pnt p2) {
        return new Pnt(p1.x - p2.x, p1.y - p2.y);
    }
    public static double norm(Pnt p) {
        return Math.sqrt(p.x * p.x + p.y * p.y);
    }
    public static int norm2(Pnt p) {
        return (p.x * p.x + p.y * p.y);
    }
    public static int dot(Pnt p1, Pnt p2) {
        return p1.x * p2.x + p1.y * p2.y;
    }
    public static int cross(Pnt p1, Pnt p2) {
        return p1.x * p2.y - p1.y * p2.x;
    }
    public static double dist(Pnt p1, Pnt p2) {
        return norm(substr(p1, p2));
    }
    public static int dist2(Pnt p1, Pnt p2) {
        return norm2(substr(p1, p2));
    }
}

// ------------- class Edge ------------------------------
class Edge {
    public Pnt p1, p2, vect;    //vector p1 -> p2
    public double norm;
    public Edge() {};
    public Edge(Pnt p1n, Pnt p2n) {
        p1 = p1n;
        p2 = p2n;
        vect = G2D.substr(p2, p1);
        norm = G2D.norm(vect);
    }
    public Edge(int x1, int y1, int x2, int y2) {
        p1 = new Pnt(x1, y1);
        p2 = new Pnt(x2, y2);
        vect = G2D.substr(p2, p1);
        norm = G2D.norm(vect);
    }
    boolean eq(double a, double b){
        return Math.abs(a-b) < 1e-9;
    }
    // ---------------------------------------------------
    public boolean intersect(Edge other) {
        //do edges "this" and "other" intersect?
        if (Math.min(p1.x,p2.x) > Math.max(other.p1.x,other.p2.x)) return false;
        if (Math.max(p1.x,p2.x) < Math.min(other.p1.x,other.p2.x)) return false;
        if (Math.min(p1.y,p2.y) > Math.max(other.p1.y,other.p2.y)) return false;
        if (Math.max(p1.y,p2.y) < Math.min(other.p1.y,other.p2.y)) return false;

        int den = other.vect.y*vect.x-other.vect.x*vect.y;
        int num1 = other.vect.x*(p1.y-other.p1.y)-other.vect.y*(p1.x-other.p1.x);
        int num2 =       vect.x*(p1.y-other.p1.y)-      vect.y*(p1.x-other.p1.x);

        //parallel edges
        if (den==0)
        {   if (Math.min(other.dist2(this),dist2(other)) > 0)
            return false;
            //on the same line - "not intersect" only if one of the vertices is common,
            //and the other doesn't belong to the line
            if ((this.p1==other.p1 && eq(G2D.dist(this.p2, other.p2) , this.norm + other.norm)) ||
                    (this.p1==other.p2 && eq(G2D.dist(this.p2, other.p1) , this.norm + other.norm)) ||
                    (this.p2==other.p1 && eq(G2D.dist(this.p1, other.p2) , this.norm + other.norm)) ||
                    (this.p2==other.p2 && eq(G2D.dist(this.p1, other.p1) , this.norm + other.norm)))
                return false;
            return true;
        }

        //common vertices
        if (this.p1==other.p1 || this.p1==other.p2 || this.p2==other.p1 || this.p2==other.p2)
            return false;

        double u1 = num1*1./den;
        double u2 = num2*1./den;
        if (u1<0 || u1>1 || u2<0 || u2>1)
            return false;
        return true;
    }
    // ---------------------------------------------------
    public double dist(Pnt p) {
        //distance from p to the edge
        if (G2D.dot(vect, G2D.substr(p, p1))<=0)
            return G2D.dist(p, p1);         //from p to p1
        if (G2D.dot(vect, G2D.substr(p, p2))>=0)
            return G2D.dist(p, p2);         //from p to p2
        //distance to the line itself
        return Math.abs(-vect.y*p.x+vect.x*p.y+p1.x*p2.y-p1.y*p2.x)/norm;
    }
    // ---------------------------------------------------
    public double dist2(Edge other) {
        //distance from the closest of the endpoints of "other" to "this"
        return Math.min(dist(other.p1), dist(other.p2));
    }
}

// ------------- class SmallPolygon itself --------------
class SmallPolygonsVis {
    final int SZ = 700;             // field size
    int NP, N, Npoly;               // number of points given, max number of polygons and number of polygons selected
    Pnt[] p;                        // coordinates of points (fixed)
    int[] pointsPar;                // coordinates of points (as an array parameter)
    int[][] polys;                  // indices of points which form polygons
    int[] polysVert;                // number of vertices in each poly
    boolean valid[];
    int[] used;                     // which poly uses this point?
    HashSet<Integer> badEdges = new HashSet<>(); // intersecting edges
    // ---------------------------------------------------
    void generate(String seed) {
        try {
            SecureRandom rnd = SecureRandom.getInstance("SHA1PRNG");
            rnd.setSeed(Long.parseLong(seed));
            // generate points by sampling each coordinate uniformly, without duplicates
            int i, j, k;
            // number of points
            if (seed.equals("1"))
                NP = 10;
            else {
                int testSize = rnd.nextInt(3);
                if (testSize == 0) NP = rnd.nextInt(80) + 20;
                else if (testSize == 1) NP = rnd.nextInt(400) + 100;
                else NP = rnd.nextInt(1001) + 500;
            }
            System.out.println("NP = " + NP);
            p = new Pnt[NP];

            // generate the points
            boolean ok;
            for (i = 0; i < NP; ++i)
            {
                do {
                    p[i] = new Pnt(rnd.nextInt(SZ), rnd.nextInt(SZ));
                    ok = true;
                    for (j = 0; j < i && ok; ++j)
                        if (p[i].equals(p[j]))
                            ok = false;
                }
                while (!ok);
            }

            // convert points to parameter array
            pointsPar = new int[2 * NP];
            for (i = 0; i < NP; ++i) {
                pointsPar[2 * i] = p[i].x;
                pointsPar[2 * i + 1] = p[i].y;
            }

            if (manual) {
                // and to coordToPoint
                coordToPoint = new int[SZ][SZ];
                for (i = 0; i < SZ; ++i)
                    Arrays.fill(coordToPoint[i], -1);
                int x, y;
                for (i = 0; i < NP; ++i)
                    for (j = -1; j <= 1; ++j)
                        for (k = -1; k <= 1; ++k)
                        {
                            x = p[i].x + j;
                            y = p[i].y + k;
                            if (x >= 0 && x < SZ && y >= 0 && y < SZ)
                                coordToPoint[x][y] = i;
                        }
            }

            N = rnd.nextInt(19) + 2;
            if (seed.equals("1"))
                N = 3;
            System.out.println("N = " + N);
        }
        catch (Exception e) {
            addFatalError("An exception occurred while generating test case.");
            e.printStackTrace();
        }
    }
    // ---------------------------------------------------
    String validatePoly(int[] poly, int n) {
        // check that the polygon satisfies all individual conditions
        if (n < 3)
            return "a polygon must have at least 3 vertices.";

        // simple polygon: no self-intersections except in common vertices of edges
        for (int i = 0; i < n; ++i)
            for (int j = i + 1; j < n; ++j) {
                // check intersection of i..i+1 and j..j+1
                Edge e1 = new Edge(p[poly[i]], p[poly[(i + 1) % n]]);
                Edge e2 = new Edge(p[poly[j]], p[poly[(j + 1) % n]]);
                if (e1.intersect(e2)) {
                    badEdges.add(poly[i]);
                    badEdges.add(poly[j]);
                    return "edges "+poly[i]+"-"+poly[(i+1)%n]+" and "+poly[j]+"-"+poly[(j+1)%n]+" intersect";
                }
            }
        return "";
    }
    // ---------------------------------------------------
    double area(int[] poly, int n) {
        // trapezium method
        double s = 0;
        for (int i = 0; i < n; i++)
            s += (p[poly[(i + 1) % n]].y + p[poly[i]].y) * (p[poly[(i + 1) % n]].x - p[poly[i]].x) / 2.0;
        return Math.abs(s);
    }
    // ---------------------------------------------------
    double calcScore() {
        // calculate the score of current set of polygons (sum of areas), including full validity check
        // will be called from interactive editing to show the results of changes
        // 1. there are at most N polygons
        if (Npoly > N) {
            addFatalError("You can have at most " + N + " polygons.");
            return 0;
        }

        // 2. each point is used by one of polygons (no polygons using same point checked earlier)
        for (int i = 0; i < used.length; ++i)
            if (used[i] == -2)
            {
                addFatalError("Point " + i + " is not used in any polygon.");
                return 0;
            }

        // 3. each polygon is valid on its own
        for (int i = 0; i < polys.length; ++i)
        {
            if (manual && polysVert[i] == 0)
            {
                // polysVert[i] = 0 in manual mode means this is a deleted polygon; check only non-deleted ones
                continue;
            }
            if (!valid[i] && strict)
            {
                addFatalError("Polygon " + i + " is not valid: " + validatePoly(polys[i], polysVert[i]));
                return 0;
            }
        }

        // 4. no two polygons intersect
        // for an intersection to be detected, polygons have to have at least 2 vertices, so deleted polygons in manual mode have no effect
        for (int i = 0; i < polys.length; ++i)
            for (int j = 0; j < polysVert[i]; ++j) {
                for (int k = i + 1; k < polys.length; ++k)
                    for (int l = 0; l < polysVert[k]; ++l) {
                        // check intersection of edge j..j+1 of polygon i and edge l..l+1 of polygon k
                        Edge e1 = new Edge(p[polys[i][j]], p[polys[i][(j + 1) % polysVert[i]]]);
                        Edge e2 = new Edge(p[polys[k][l]], p[polys[k][(l + 1) % polysVert[k]]]);
                        if (e1.intersect(e2)) {
                            badEdges.add(polys[i][j]);
                            badEdges.add(polys[k][l]);
                            addFatalError("edges "+polys[i][j]+"-"+polys[i][(j + 1) % polysVert[i]]+" and "+polys[k][l]+"-"+polys[k][(l + 1) % polysVert[k]]+" intersect");
                            return 0;
                        }
                    }
            }

        // now, if all are valid, score is always non-0
        double score = 0;
        for (int i = 0; i < polys.length; ++i)
            score += area(polys[i], polysVert[i]);
        return score;
    }
    // ------------- visualization part ----------------------
    static String exec;
    static boolean vis, manual, debug, strict;
    static Process proc;
    JFrame jf;
    Vis v;
    InputStream is;
    OutputStream os;
    BufferedReader br;
    // problem-specific drawing params
    final int SZX = SZ+2+100,SZY=SZ+2;
    volatile boolean ready;
    volatile int Ncur;
    volatile int[] Pcur;
    int[][] coordToPoint;
    // ---------------------------------------------------
    String[] choosePolygons(int[] points, int N) throws IOException
    {   // pass the params to the solution and get the return
        int i;
        StringBuffer sb = new StringBuffer();
        sb.append(points.length).append('\n');
        for (i = 0; i < points.length; ++i)
            sb.append(points[i]).append('\n');
        sb.append(N).append('\n');
        os.write(sb.toString().getBytes());
        os.flush();
        // get the return - an array of strings
        int nret = Integer.parseInt(br.readLine());
        System.out.println(nret);
        String[] ret = new String[nret];
        for (i = 0; i < nret; ++i)
            ret[i] = br.readLine();
        return ret;
    }
    // ---------------------------------------------------
    void draw() {
        if (!vis) return;
        v.repaint();
    }
    // ---------------------------------------------------
    public class Vis extends JPanel implements MouseListener, WindowListener {
        public void paint(Graphics g) {
            try {
                //do painting here
                int i,j,n;
                char[] ch;
                BufferedImage bi = new BufferedImage(SZX+10,SZY+10,BufferedImage.TYPE_INT_RGB);
                Graphics2D g2 = (Graphics2D)bi.getGraphics();
                //background
                g2.setColor(new Color(0xD3D3D3));
                g2.fillRect(0,0,SZX+10,SZY+10);
                g2.setColor(Color.WHITE);
                g2.fillRect(0,0,SZ+1,SZ+1);
                //frame
                g2.setColor(Color.BLACK);
                g2.drawRect(0,0,SZ+1,SZ+1);

                //sides
                for (i=0; i<polys.length; i++) {
                    n = polysVert[i];
                    if (valid[i]) {
                        float hue = (float)(i) / polys.length;
                        g2.setColor(Color.getHSBColor(hue, 0.9f, 1.0f));
                        int[] xPoints = new int[n];
                        int[] yPoints = new int[n];
                        for (j=0; j<n; j++) {
                            xPoints[j] = p[polys[i][j]].x;
                            yPoints[j] = (SZ-1-p[polys[i][j]].y);
                        }
                        g2.fillPolygon(xPoints, yPoints, n);
                    }
                    if (valid[i])
                        g2.setColor(Color.GREEN);
                    else
                        g2.setColor(Color.RED);
                    for (j=0; j<n; j++) {
                        g2.drawLine(p[polys[i][j]].x, (SZ-1-p[polys[i][j]].y), p[polys[i][(j+1)%n]].x, (SZ-1-p[polys[i][(j+1)%n]].y));
                    }
                    if (badEdges.size()>0) {
                        g2.setColor(Color.RED);
                        g2.setStroke(new BasicStroke(3));
                        for (j=0; j<n; j++) {
                            if (badEdges.contains(polys[i][j]))
                                g2.drawLine(p[polys[i][j]].x, (SZ-1-p[polys[i][j]].y), p[polys[i][(j+1)%n]].x, (SZ-1-p[polys[i][(j+1)%n]].y));
                        }
                        g2.setStroke(new BasicStroke(1));
                    }
                }
                //draw current poly
                g2.setColor(new Color(0x6495ED));
                for (i=0; i<Ncur; i++) {
                    g2.drawLine(p[Pcur[i]].x, (SZ-1-p[Pcur[i]].y), p[Pcur[(i+1)%Ncur]].x, (SZ-1-p[Pcur[(i+1)%Ncur]].y));
                }

                //"buttons"
                if (manual) {
                    g2.setColor(Color.BLACK);
                    ch = ("SUBMIT").toCharArray();
                    g2.setFont(new Font("Arial",Font.BOLD,16));
                    g2.drawChars(ch,0,ch.length,SZ+20,30);
                    g2.drawRect(SZ+12,8,90,30);

                    ch = ("ADD POLY").toCharArray();
                    g2.setFont(new Font("Arial",Font.BOLD,14));
                    g2.drawChars(ch,0,ch.length,SZ+18,109);
                    g2.drawRect(SZ+12,88,90,30);

                    ch = ("DEL POLY").toCharArray();
                    g2.setFont(new Font("Arial",Font.BOLD,14));
                    g2.drawChars(ch,0,ch.length,SZ+19,149);
                    g2.drawRect(SZ+12,128,90,30);
                }

                //current score
                ch = (""+calcScore()).toCharArray();
                g2.setFont(new Font("Arial",Font.BOLD,14));
                g2.drawChars(ch,0,ch.length,SZ+10,200);

                //points with small digits near them
                g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
                g2.setFont(new Font("Arial",Font.PLAIN,10));
                for (i=0; i<NP; i++)
                {   if (used[i]>-1) {
                    if (valid[used[i]])
                        g2.setColor(Color.GREEN);
                    else g2.setColor(Color.RED);
                }
                else {
                    if (used[i]==-1) {
                        //a special highlight for last point in the polygon
                        if (Pcur[Ncur-1]==i)
                            g2.setColor(new Color(0x6495ED));
                        else g2.setColor(new Color(0x000080));
                    }
                    else g2.setColor(Color.BLACK);
                }
                    g2.fillOval(p[i].x-2, SZ-1-p[i].y-2,5,5);
                    if (debug) {
                        g2.setColor(Color.BLACK);
                        ch = (i+"").toCharArray();
                        g2.drawChars(ch,0,ch.length, p[i].x+2, SZ-1-p[i].y-2);
                    }
                }

                g.drawImage(bi,0,0,SZX+10,SZY+10,null);
            }
            catch (Exception e) { e.printStackTrace(); }
        }
        public Vis() {
            if (manual)
                addMouseListener(this);
            jf.addWindowListener(this);
        }
        // ---------------------------------------------------
        //MouseListener
        public void mouseClicked(MouseEvent e) {
            //identify the buttons or the click on the field
            int x = e.getX(), y = e.getY(), i, j;
            if (x>SZ) {
                //can be only mode modifiers
                if (y >= 8 && y <= 38) {
                    //"SUBMIT"
                    ready = true;
                }
                if (y >= 88 && y <= 118) {
                    //"ADD POLY"
                    String val = validatePoly(Pcur, Ncur);
                    if (val.length()!=0) {
                        System.out.println("Current polygon is invalid: "+val);
                    }
                    else {
                        //actually commit the polygon - find first slot with no poly in it and add it there
                        for (i=0; i<polys.length; i++)
                            if (polysVert[i]==0)
                                break;
                        if (i == polys.length) {
                            System.out.println("Can't have more than "+polys.length+" polygons.");
                        } else {
                            //put current poly to the slot
                            if (debug) System.out.println("Adding current polygon to slot "+i);
                            polysVert[i] = Ncur;
                            valid[i] = true;    //already verified
                            if (polys[i] == null || polys[i].length < Ncur)
                                polys[i] = new int[Ncur];
                            for (j=0; j<Ncur; j++) {
                                polys[i][j] = Pcur[j];
                                used[Pcur[j]] = i;
                            }
                            Ncur = 0;
                            Npoly++;
                        }
                    }
                }
                if (y >= 128 && y <= 158) {
                    //"DEL POLY"
                    //delete the currently selected poly (and unmark used points)
                    if (debug) System.out.println("Deleting current polygon");
                    for (i=0; i<Ncur; i++)
                        used[Pcur[i]] = -2;
                    Ncur = 0;
                }
                draw();
                return;
            }
            //now, the clicks weren't buttons - they were points' locations
            y = SZ-y-1;
            int indP = coordToPoint[x][y], indPoly;
            if (indP==-1)
                return;
            //now process the point
            indPoly = used[indP];

            //three scenarios: adding a point to current poly, removing it or choosing a poly to be edited
            if (Ncur == 0 && indPoly > -1) {
                //this polygon is moved to editing: delete it from the list and free its parameters
                if (debug) System.out.println("Editing polygon "+indPoly);
                Ncur = polysVert[indPoly];
                polysVert[indPoly] = 0;
                valid[indPoly] = false;
                for (i=0; i<Ncur; i++) {
                    Pcur[i] = polys[indPoly][i];
                    polys[indPoly][i] = -1;
                    used[Pcur[i]] = -1;
                }
                Npoly--;
            } else
            if (Ncur > 0 && indPoly == -1 && Pcur[Ncur-1]==indP) {
                //this point is last in the polygon - remove it
                if (debug) System.out.println("Removing point "+indP+" from current polygon");
                Ncur --;
                used[indP] = -2;
                Pcur[Ncur] = -2;
            } else
            if (indPoly == -2) {
                if (debug) System.out.println("Adding point "+indP+" to current polygon");
                Pcur[Ncur] = indP;
                Ncur++;
                used[indP] = -1;
            } else {
                if (debug) System.out.println("Invalid action");
                return;
            }

            draw();
        }
        public void mousePressed(MouseEvent e) { }
        public void mouseReleased(MouseEvent e) { }
        public void mouseEntered(MouseEvent e) { }
        public void mouseExited(MouseEvent e) { }
        //WindowListener
        public void windowClosing(WindowEvent e) {
            if(proc != null)
                try { proc.destroy(); }
                catch (Exception ex) { ex.printStackTrace(); }
            System.exit(0);
        }
        public void windowActivated(WindowEvent e) { }
        public void windowDeactivated(WindowEvent e) { }
        public void windowOpened(WindowEvent e) { }
        public void windowClosed(WindowEvent e) { }
        public void windowIconified(WindowEvent e) { }
        public void windowDeiconified(WindowEvent e) { }
    }

    public SmallPolygonsVis() {}

    void addFatalError(String message) {
        System.out.println(message);
    }
}

