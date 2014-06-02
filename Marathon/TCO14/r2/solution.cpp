#include <iostream>
#include <sstream>
#include <cstdio>
#include <vector>
#include <string>
#include <algorithm>
#include <set>
#include <map>
#include <cmath>
#include <unordered_set>
#include <unordered_map>
#include <queue>
#include <utility>
#include <cassert>
#include <numeric>
#include <sstream>
#include <limits>
using namespace std;

#define REQUIRE(cond, message) \
    do { \
        if (!(cond)) { \
            std::cerr << message << std::endl; \
            assert(false); \
        } \
    } while (false)

#define forn(i, n) for (int i = 0; i < int(n); ++i)
#define for1(i, n) for (int i = 1; i <= int(n); ++i)
#define forv(i, v) forn(i, (v).size())
#define pb push_back
#define mp make_pair
#define all(v) (v).begin(), (v).end()

typedef vector<int> vi;
typedef vector<vi> vvi;
typedef long long ll;
typedef vector<ll> vl;
typedef pair<int, int> pii;
typedef vector<string> vs;
typedef long double ld;


const int NMAX_COORD = 1000000;

enum Side { TOP, LEFT, BOTTOM, RIGHT };

class RectanglesAndHoles
{
public:
    vi place(vi a, vi b);
};

struct Block
{
    int height;
    int width;
    int id;
};

struct Place
{
    int x;
    int y;
    int orient;
    int id;
};

struct Segment
{
    int length;
    int id;
};

typedef vector<Segment> Segments;

vector<Block> toBlocks(const vi& a, const vi& b)
{
    vector<Block> blocks(a.size());
    forv(i, a) {
        Block& block = blocks[i];
        block.width = a[i];
        block.height = b[i];
        block.id = i;
    }
    return blocks;
}

vi toResultFormat(const vector<Place>& placement)
{
    vi result(placement.size() * 3);
    forv(i, placement) {
        const Place& p = placement[i];
        result[3 * i] = p.x;
        result[3 * i + 1] = p.y;
        result[3 * i + 2] = p.orient;
        assert(0 <= p.orient && p.orient <= 1);
    }
    return result;
}

int calcMaxLength(const Segments& segments)
{
    int length = 0;
    for (const Segment& segment : segments) {
        length += segment.length;
    }
    return length;
}

void place(Side side, int x, int y,
           const Segments& segments,
           const vector<Block>& blocks,
           vector<Place>* placement)
{
    for (const Segment& segment : segments) {
        int index = segment.id;
        bool isHorizontal = side == Side::TOP || side == Side::BOTTOM;
        const Block& block = blocks[index];

        int orient = (isHorizontal ? block.width : block.height) != segment.length;
        int xShift = 0;
        int yShift = 0;
        if (side == Side::BOTTOM) {
            yShift -= orient ? block.width : block.height;
        }
        else if (side == Side::LEFT) {
            xShift -= orient ? block.height : block.width;
        }

        Place place = {x + xShift, y + yShift, orient, block.id};
        placement->push_back(place);
        if (isHorizontal) {
            x += segment.length;
        }
        else {
            y += segment.length;
        }
    }
}


void buildRectangleGreedy(const Segments& segments,
                          Segments* top,
                          Segments* left,
                          Segments* bottom,
                          Segments* right)
{
    vector<Segments*> choices {top, left, bottom, right};
    vi len(choices.size());
    Segments sorted = segments;
    sort(all(sorted), [](const Segment& a, const Segment& b) { return a.length > b.length; });
    for (const Segment& segment : sorted) {
        int minLenId = 0;
        for (size_t i = 1; i < len.size(); ++i) {
            if (len[i] < len[minLenId]) minLenId = i;
        }
        len[minLenId] += segment.length;
        choices[minLenId]->push_back(segment);
    }
}
                          

void buildMaxAreaRectangle(const Segments& segments,
                           Segments* top,
                           Segments* left,
                           Segments* bottom,
                           Segments* right)
{
    buildRectangleGreedy(segments, top, left, bottom, right);
}

vector<Place> placeMaxArea(const vector<Block>& blocks)
{
    int numBlocks = blocks.size();
    Segments segments(numBlocks);
    forn(i, numBlocks) {
        segments[i].id = i;
        segments[i].length = max(blocks[i].width, blocks[i].height);
    }
    Segments topSegments, bottomSegments, leftSegments, rightSegments;
    buildMaxAreaRectangle(segments, &topSegments, &leftSegments,
                          &bottomSegments, &rightSegments);
    int topLength = calcMaxLength(topSegments);
    int bottomLength = calcMaxLength(bottomSegments);
    int leftLength = calcMaxLength(leftSegments);
    int rightLength = calcMaxLength(rightSegments);
    if (bottomLength < topLength) {
        swap(bottomLength, topLength);
        bottomSegments.swap(topSegments);
    }
    if (leftLength < rightLength) {
        swap(leftLength, rightLength);
        leftSegments.swap(rightSegments);
    }
    vector<Place> placement;
    place(Side::LEFT, 0, 0, leftSegments, blocks, &placement);
    place(Side::RIGHT, topLength, 0, rightSegments, blocks, &placement);
    place(Side::BOTTOM, 0, 0, bottomSegments, blocks, &placement);
    place(Side::TOP, 0, rightLength, topSegments, blocks, &placement);
    return placement;
}

vi RectanglesAndHoles::place(vi a, vi b)
{
    vector<Block> blocks = toBlocks(a, b);
    vector<Place> placement = placeMaxArea(blocks);
    return toResultFormat(placement);
}

int main()
{
    ios_base::sync_with_stdio(false);
    RectanglesAndHoles task;
    int n;
    cin >> n;
    vi a(n);
    forn(i, n) cin >> a[i];
    vi b(n);
    forn(i, n) cin >> b[i];
    vi result = task.place(a, b);
    forv(i, result) {
        cout << result[i] << endl;
    }
    return 0;
}
