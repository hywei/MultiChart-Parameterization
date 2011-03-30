/* ================== Library Information ================== */
// [Name] 
// MeshLib Library
//
// [Developer]
// Xu Dong
// State Key Lab of CAD&CG, Zhejiang University
// 
// [Date]
// 2005-08-05
//
// [Goal]
// A general, flexible, versatile and easy-to-use mesh library for research purpose.
// Supporting arbitrary polygonal meshes as input, but with a 
// primary focus on triangle meshes.

/* ================== File Information ================== */
// [Name]
// Macro.h
//
// [Developer]
// Xu Dong
// State Key Lab of CAD&CG, Zhejiang University
// 
// [Date]
// 2005-08-05
//
// [Goal]
// Defining the global macros and enumerations.



#pragma once


#define PLATFORM_WINDOWS 1

#ifdef PLATFORM_WINDOWS
#define MESHLIB_ASSERT(a) ASSERT(a)
#endif

#include <cmath>

/* ================== Enumerations ================== */

typedef enum {SOLID_SMOOTH, SOLID_FLAT, WIREFRAME, SOLID_AND_WIREFRAME, VERTICES, HIGHLIGHT_ONLY, TEXTURE_MAPPING,
	SELECT_NULL, SELECT_VERTEX, SELECT_EDGE, SELECT_FACE, } RENDER_MODE;
 
typedef enum {TRANSLATE, SCALE, ROTATE, NONE} TRANSFORM_MODE;

typedef enum {FACE, EDGE, VERTEX} PICKING_MODE;

typedef enum {AXIS_X, AXIS_Y, AXIS_Z} AXIS;

typedef enum {R=0, G, B} COLOR;



/* ================== Constants ================== */

const double PI = 3.1415926535897932384626433832795;
const double PI_12 = 1.5707963267948966192313216916395;
const double LARGE_ZERO_EPSILON = 1.0E-8;
const double SMALL_ZERO_EPSILON = 1.0E-16;



/* ================== Math Macros ================== */

// Min & Max & ABS
template <typename T> inline T Max3(const T& a, const T& b, const T& c)
{
	return (a > b) ? ( a > c ? a : c) : ( b > c ? b : c);
}
template <typename T> inline T Min3(const T& a, const T& b, const T& c)
{
	return (a < b) ? ( a < c ? a : c) : ( b < c ? b : c);
}


// Equal within epsilon tolerance
inline bool ALMOST_EQUAL_LARGE(double x, double y)
{
	return (fabs(x - y) < LARGE_ZERO_EPSILON);
}
inline bool ALMOST_EQUAL_SMALL(double x, double y)
{
	return (fabs(x - y) < SMALL_ZERO_EPSILON);
}

inline bool GreaterEqual(double a, double b, double epsilon=LARGE_ZERO_EPSILON)
{
	return (a > b) || ( fabs(a - b) < epsilon);
}
inline bool LessEqual(double a, double b, double epsilon = LARGE_ZERO_EPSILON)
{
	return (a < b) || ( fabs(a - b) < epsilon);
}

// Degree & Radius
#define DEG2RAD(d) ((d)*0.017453292519943295769236907684886)
#define RAD2DEG(r) ((r)*57.29577951308232087679815481410500)

