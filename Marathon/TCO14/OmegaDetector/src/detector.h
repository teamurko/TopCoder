#include <string.h>
#include <ctype.h>
#include <math.h>

// define the maximum size image
#define FULL_IMAGE_MAX_COLS 1280
#define FULL_IMAGE_MAX_ROWS 960
//define the maximum result length
#define MAX_RESULT_LENGTH   4096
// define maximum number of barcode error corrections
#define MAX_EC_CORRECTIONS  4096
// define maximum number of regions-of-interest
#define MAX_ROI_COUNT       6

// Types - Must use the following types to properly target
// the desired machine architecture:
// U8 - means an unsigned 8 bit variable
// S8 - means a signed 8 bit variable
// U16 - means an unsigned 16 bit variable
// S16 - means an signed 16 bit variable
// U32 - means an unsigned 32 bit variable
// S32 - means an signed 32 bit variable
// U64 - means an unsigned 64 bit variable
// S64 - means an signed 64 bit variable
// F32 - means a single precision float
// F64 - means a double precision float

// All variable names should be prefixed in such a way as
// to indicate the type of that variable. This is sometimes
// referred to as Hungarian notation. We have heard all of the
// arguments for and against Hungarian notation. Don't fuss or
// complain about it - just do it.

// On many systems the following typedefs can be used to
// specify variable types of the correct size:
typedef unsigned char  U8;        // variable name prefix is: u8
typedef char   S8;                // variable name prefix is: s8
typedef unsigned short U16;       // variable name prefix is: u16
typedef short  S16;               // variable name prefix is: s16
typedef unsigned long  U32;       // variable name prefix is: u32
typedef long   S32;               // variable name prefix is: s32
typedef unsigned long long U64;   // variable name prefix is: u64
typedef long long S64;            // variable name prefix is: s64
typedef float  F32;               // variable name prefix is: f32
typedef double F64;               // variable name prefix is: f64

// All pixel manipulation functions should use
// variable type PIXEL_T and no other built-in type
// or typedef'ed type.
typedef unsigned char PIXEL_T;   // variable name prefix is: px

// Use this structure to return a point
// characterized as X,Y coordinates within
// an image.
typedef struct {
U16 u16X;             // X-coordinate of the point (also the column number)
U16 u16Y;             // Y-coordinate of the point (also the row number)
} st_POINT;                       // variable name prefix is: st

// An image is stored in memory as a 2D
// array, of type PIXEL_T, in row
// major order.  In addition to rows
// and columns, stride is specified as the
// the distance in memory between two rows
// of pixels at the same column location.
// This allows the easy identification of a
// subregion of a larger image.  Of course,
// the full image can be identified in exactly
// the same way.
typedef struct {
U16 u16Rows;          // number of rows in the image region
U16 u16Cols;          // number of columns in the image region
U16 u16Stride;        // distance in memory between rows
PIXEL_T * pxImage; // pointer to the first pixel of the region
} st_IMAGE_REGION;             // variable name prefix is: st

// The following definitions are used to indicate which barcode symbology
// or multiple symbologies are present in an image region.
enum { en1DDetect=0,en2DDetect };
#define SYMBOLOGY_MASK(s) (1 << (s))

// Results of a call to Topcoder_ROI_Finder function are
// returned in an array of the following type. The array
// of st_ROI_INFO structures exists before
// Topcoder_ROI_Finder functions is called, so the function
// need only populate the structure fields for each actual result.
// and zero the scores for each non result.
typedef struct {
S32 s32ConfidenceScore; // score proportional to how much the region looks like a barcode
st_POINT stMinXY;     // minimum (X,Y) (column, row) coordinates of the rectangular region of interest
st_POINT stMaxXY;     // maximum (X,Y) (column, row) coordinates of the rectangular region of interest
U32 u32Symbologies;   // region may be tagged to indicate barcode symbologies
} st_ROI_INFO;


// You can add your code here .....

S8 Topcoder_ROI_Finder(st_IMAGE_REGION *pstImageRegion, st_ROI_INFO *pstDetections)
{
    // You can also add your code here .....
    // populate pstDetections and return the number of barcodes detected.

    return 0;
}




class OmegaDetector
{
public:
    std::vector<int> ROI_Finder(const std::vector<int> &imageData, int numOfRows, int numOfCols)
    {
        std::vector<int> ROIs;

        S8 s8RegionCount=0;
        st_ROI_INFO stInterestingRegions[MAX_ROI_COUNT];
        st_IMAGE_REGION stRegion;

        stRegion.u16Rows = numOfRows;
        stRegion.u16Cols = numOfCols;
        stRegion.u16Stride = numOfCols;
        stRegion.pxImage = (PIXEL_T*)(&imageData[0]);

        s8RegionCount = Topcoder_ROI_Finder(&stRegion, &stInterestingRegions[0]);

        ROIs.resize(6*s8RegionCount);
        int j = 0;
        for (int i=0;i<s8RegionCount;i++)
        {
            ROIs[j++] = stInterestingRegions[i].stMinXY.u16X;
            ROIs[j++] = stInterestingRegions[i].stMinXY.u16Y;
            ROIs[j++] = stInterestingRegions[i].stMaxXY.u16X;
            ROIs[j++] = stInterestingRegions[i].stMaxXY.u16Y;
            ROIs[j++] = stInterestingRegions[i].u32Symbologies;
            ROIs[j++] = stInterestingRegions[i].s32ConfidenceScore;
        }

        return ROIs;
    }
};
