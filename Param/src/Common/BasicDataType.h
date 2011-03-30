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
// BasicDataType.h
//
// [Developer]
// Xu Dong
// State Key Lab of CAD&CG, Zhejiang University
// 
// [Date]
// 2005-08-05
//
// [Goal]
// Defining the basic data type used in the MeshLib library.
// Including coordinate, normal, color, texture coordinate, various indices, etc

#pragma once

#include "Macro.h"
#include <cmath>
#include <vector>
#include <string>
#include <cstdio>


/* ================== Various Indices ================== */
typedef int VertexID;    // Vertex index
typedef int	EdgeID;      // HalfEdge index
typedef int	FaceID;      // Face index

typedef int	PointID;     // Point index
typedef int	LineID;      // Line index
typedef int	PolygonID;   // Polygon index

typedef int BdyID;       // Boundary index
typedef int ComponentID; // Component index

typedef int	ModelID;     // Model index
typedef int	SceneID;     // Scene index



/* ================== Coordinate (3D) ================== */
class Coord
{ 
private:
    double pos[3];

public:
    // Constructor
    Coord() 
    { 
        pos[0] = pos[1] = pos[2] = 0.0; 
    }
	Coord(double v[3])
    {
        pos[0] = v[0];
        pos[1] = v[1];
        pos[2] = v[2];
    }

	Coord(const Coord &c)
    {
        pos[0] = c.pos[0];
        pos[1] = c.pos[1];
        pos[2] = c.pos[2];
    }

	Coord(double _x, double _y = 0.0, double _z = 0.0)
    {
        pos[0] = _x;
        pos[1] = _y;
        pos[2] = _z;
    }

    // Destructor
    ~Coord()
    {
        
    }

    // Initializer
    inline void setCoords(double _x = 0.0, double _y = 0.0, double _z = 0.0);

    // Basic operator
	inline Coord operator +(const Coord &) const;
	inline Coord operator -(const Coord &) const;
	inline Coord operator *(double) const;
	inline Coord operator /(double) const;
	inline Coord operator -() const;
    
    inline Coord& operator +=(const Coord &);
	inline Coord& operator -=(const Coord &);
	inline Coord& operator *=(double);
	inline Coord& operator /=(double);

    inline Coord& operator =(const Coord &);

	inline bool operator ==(const Coord &) const;
	inline bool operator !=(const Coord &) const;

	inline double& operator [](size_t) ;
	inline double operator [](size_t) const;

    // Advanced operator
	inline friend double dot(const Coord&, const Coord&);    
    inline friend double angle(const Coord&, const Coord&);
	inline friend Coord cross(const Coord&, const Coord&);
    
	inline void print(std::string name);

	inline double abs() const;
	inline double sqrabs() const;

	inline Coord unit() const;
    inline bool normalize();

    // Transformer
	inline Coord Spin(Coord SpinAxes, Coord Center, double SpinAngle);
	inline Coord Scale(Coord scale);
	inline Coord Rotate(double delta_theta, double delta_phi, double theta);

	//
	inline double x() const;
	inline double y() const;
	inline double z() const;

	//
};



/* ================== Coordinate (2D) ================== */
class Coord2D 
{ 
private:
    double pos[2];

public:
    // Constructor
    Coord2D()
    {
        pos[0] = pos[1] = 0.0;
    }
    
    Coord2D(double _x, double _y)
    {
        pos[0] = _x;
        pos[1] = _y;
    }
    
    Coord2D(double v[2])
    {
        pos[0] = v[0];
        pos[1] = v[1];
    }
    
    Coord2D(const Coord2D &c)
    {
        pos[0] = c.pos[0];
        pos[1] = c.pos[1];
    }
    
    // Destructor
    ~Coord2D()
    {
        
    }

    // Initializer
	inline void setCoords(double _x = 0.0, double _y = 0.0);

    // Basic operator
	inline Coord2D operator +(const Coord2D &) const;
	inline Coord2D operator -(const Coord2D &) const;
	inline Coord2D operator *(double) const;
	inline Coord2D operator /(double) const;
	inline Coord2D operator -() const;

	inline Coord2D& operator +=(const Coord2D &);
	inline Coord2D& operator -=(const Coord2D &);
	inline Coord2D& operator *=(double);
	inline Coord2D& operator /=(double);

    inline Coord2D& operator =(const Coord2D &);

	inline bool operator ==(const Coord2D &) const;
	inline bool operator !=(const Coord2D &) const;

	inline double& operator [](int);
	inline double operator [](int) const;

    // Advanced operator
	inline friend double dot(const Coord2D&, const Coord2D&);
	inline friend double angle(const Coord2D&, const Coord2D&);

	inline double abs() const;
	inline double sqrabs() const;
    
	inline Coord2D unit() const;
    inline bool normalize();

	inline double x() const;
	inline double y() const;
};



/* ================== Color (RGB) ================== */
class Color
{
private:
	double rgb[3];
    
public:
    // Constructor
    Color()
    {
        rgb[0] = rgb[1] = rgb[2] = 0.0;
    }
    
    Color(int r, int g, int b)               // r,g,b - [0,   255]
    {
        setColor(r, g, b);
    }
    
    Color(double r, double g, double b)      // r,g,b - [0.0, 1.0]
    {
        setColor(r, g, b);
    }
    
    // Destructor
    ~Color()
    {       
    }

    // Initializer
    inline void setColor(int r, int g, int b);             // r,g,b - [0,   255]
    inline void setColor(double r, double g, double b);    // r,g,b - [0.0, 1.0]
    inline void setRandomColor();
	
    //inline COLORREF Color2Ref();

    // Basic operator
	inline Color operator +(const Color &) const;
	inline Color operator -(const Color &) const;
	inline Color operator *(double) const;
	inline Color operator /(double) const;
	inline Color operator -() const;

	inline Color& operator +=(const Color &);
	inline Color& operator -=(const Color &);
	inline Color& operator *=(double);
	inline Color& operator /=(double);
    
    inline Color& operator =(const Color &color);

	inline bool operator ==(const Color &) const;
	inline bool operator !=(const Color & c) const;
    
	inline double& operator[](int i)
	{
		return rgb[i];
	}
	inline double operator[](int i) const
	{
		return rgb[i];
	}

    // Advanced operator
    inline void Clamp();   // Make sure r,g,b - [0.0, 1.0]

	//
	double r() const {return rgb[0];}
	double g() const {return rgb[1];}
	double b() const {return rgb[2];}
};

class CCurvature
{
public:
    // Constructor
    CCurvature()
	{
		m_kmax = 0;
		m_kmin = 0;
	}

	// Destructor
	~CCurvature()
	{
	}

	inline CCurvature& operator = (const CCurvature &curv);

public:
	Coord  m_direction_kmin;
	Coord  m_direction_kmax;
	double m_kmin;
	double m_kmax;
};


/* ================== Properties ================== */
typedef Coord   Normal;
typedef Coord2D TexCoord;
typedef int     Flag;
typedef int     Index;



/* ================== Property Arrays ================== */
typedef std::vector<Coord2D> Coord2DArray;
typedef std::vector<Coord>   CoordArray;
typedef std::vector<Color>   ColorArray;

typedef std::vector<bool>    BoolArray;
typedef std::vector<int>     IntArray;
typedef std::vector<double>  DoubleArray;
typedef std::vector<float>   FloatArray;

typedef std::vector<Normal>   NormalArray;
typedef std::vector<TexCoord> TexCoordArray;
typedef std::vector<Index>    IndexArray;
typedef std::vector<Flag>     FlagArray;

typedef std::vector<IndexArray >     PolyIndexArray;
typedef std::vector<CoordArray >     PolyCoordArray;
typedef std::vector<TexCoordArray >  PolyTexCoordArray;
typedef std::vector<DoubleArray >    PolyDataArray;
typedef std::vector<CCurvature>      CurvatureArray;


/* ================== Color Macros ================== */
#define RED			Color(1.0,0.0,0.0)
#define GREEN		Color(0.0,1.0,0.0)
#define BLUE		Color(0.0,0.0,1.0)
#define	YELLOW		Color(1.0,1.0,0.0)
#define CYAN		Color(0.0,1.0,1.0)
#define MAGENTA		Color(1.0,0.0,1.0)
#define BLACK		Color(0.0,0.0,0.0)
#define	WHITE		Color(1.0,1.0,1.0)
#define ORANGE		Color(1.0,0.5,0.0)
#define DARK_RED	Color(0.5,0.0,0.0)
#define LIGHT_RED	Color(1.0,0.5,0.5)
#define DARK_GREEN	Color(0.0,0.5,0.0)
#define LIGHT_GREEN	Color(0.5,1.0,0.5)
#define DARK_BLUE	Color(0.0,0.0,0.5)
#define LIGHT_BLUE	Color(0.5,0.5,1.0)

#define LIGHT_GREY  Color(0.3,0.3,0.3)
#define GREY        Color(0.5,0.5,0.5)
#define DARK_GREY   Color(0.8,0.8,0.8)



/* ================== Coordinate Macros ================== */
#define COORD_AXIS_X Coord(1.0, 0.0, 0.0)
#define COORD_AXIS_Y Coord(0.0, 1.0, 0.0)
#define COORD_AXIS_Z Coord(0.0, 0.0, 1.0)
#define COORD_ORIGIN Coord(0.0, 0.0, 0.0)



//////////////////////////////////////////////////////////////////////////
// inline implementations

/* ================== Coordinate (3D) ================== */

// Initializer
inline void Coord::setCoords(double _x, double _y, double _z)
{
	pos[0] = _x;
	pos[1] = _y;
	pos[2] = _z;
}

// Basic operator
inline Coord Coord::operator +(const Coord & c) const
{
	Coord coord;
	coord.pos[0] = pos[0]+c.pos[0];
	coord.pos[1] = pos[1]+c.pos[1];
	coord.pos[2] = pos[2]+c.pos[2];
	return coord;
}

inline Coord Coord::operator -(const Coord & c) const
{
	Coord coord;
	coord.pos[0] = pos[0]-c.pos[0];
	coord.pos[1] = pos[1]-c.pos[1];
	coord.pos[2] = pos[2]-c.pos[2];
	return coord;
}

inline Coord Coord::operator *(double v) const
{
	Coord coord;
	coord.pos[0] = pos[0]*v;
	coord.pos[1] = pos[1]*v;
	coord.pos[2] = pos[2]*v;
	return coord;
}

inline Coord Coord::operator /(double v) const
{
	Coord coord;
	if(ALMOST_EQUAL_SMALL(v,0.0))
	{
		fprintf(stdout,"Coord Error: Divided by Zero\n");
		return coord;
	}
	coord.pos[0] = pos[0]/v;
	coord.pos[1] = pos[1]/v;
	coord.pos[2] = pos[2]/v;
	return coord;
}

inline Coord Coord::operator -() const
{
	Coord coord;
	coord.pos[0] = -1.0*pos[0];
	coord.pos[1] = -1.0*pos[1];
	coord.pos[2] = -1.0*pos[2];
	return coord;
}

//
inline Coord& Coord::operator +=(const Coord & c)
{
	pos[0] += c.pos[0];
	pos[1] += c.pos[1];
	pos[2] += c.pos[2];
	return *this;
}

inline Coord& Coord::operator -=(const Coord & c)
{
	pos[0] -= c.pos[0];
	pos[1] -= c.pos[1];
	pos[2] -= c.pos[2];
	return *this;
}

inline Coord& Coord::operator *=(double v)
{
	pos[0] *= v;
	pos[1] *= v;
	pos[2] *= v;
	return *this;
}

inline Coord& Coord::operator /=(double v)
{
	if(ALMOST_EQUAL_SMALL(v,0.0))
	{
		fprintf(stdout,"Coord Error: Divided by Zero\n");
		return *this;
	}
	pos[0] /= v;
	pos[1] /= v;
	pos[2] /= v;
	return *this;
}

//
inline Coord& Coord::operator =(const Coord & c)
{
    pos[0] = c.pos[0];
    pos[1] = c.pos[1];
    pos[2] = c.pos[2];
    return *this;
}

//
inline bool Coord::operator ==(const Coord &c) const
{
	if(pos[0]!=c.pos[0] || pos[1]!=c.pos[1] || pos[2]!=c.pos[2])
		return false;
	return true;
}

inline bool Coord::operator !=(const Coord &c) const
{
	if(pos[0]!=c.pos[0] || pos[1]!=c.pos[1] || pos[2]!=c.pos[2])
		return true;
	return false;
}

//
inline double& Coord::operator [](size_t i)
{
    return pos[i];
}

inline double Coord::operator [](size_t i) const
{
    return pos[i];
}

// Advanced operator
inline double dot(const Coord& c1, const Coord& c2)
{
	return (c1.pos[0]*c2.pos[0] + c1.pos[1]*c2.pos[1] + c1.pos[2]*c2.pos[2]);
}

inline double angle(const Coord& c1, const Coord& c2)
{
	Coord tmp_c1(c1), tmp_c2(c2);
	tmp_c1.normalize(); 
	tmp_c2.normalize();

    double dot_ = dot(tmp_c1, tmp_c2);
    if(dot_ < -1.0)
        dot_ = -1.0;
    else if(dot_ > 1.0)
        dot_ = 1.0;
    
    return acos(dot_);
}

inline Coord cross(const Coord& c1, const Coord& c2)
{
	Coord c;
	c.pos[0] = c1.pos[1]*c2.pos[2] - c2.pos[1]*c1.pos[2];
	c.pos[1] = c2.pos[0]*c1.pos[2] - c1.pos[0]*c2.pos[2];
	c.pos[2] = c1.pos[0]*c2.pos[1] - c2.pos[0]*c1.pos[1];
	return c;
}

inline void Coord::print(std::string name)
{
    printf("%s = (%f, %f, %f)\n", name.c_str(), this->pos[0], this->pos[1], this->pos[2]);
}

//
inline double Coord::abs() const
{
	return sqrt(pos[0]*pos[0] + pos[1]*pos[1] + pos[2]*pos[2]);
}

inline double Coord::sqrabs() const
{
	return (pos[0]*pos[0] + pos[1]*pos[1] + pos[2]*pos[2]);
}

//
inline Coord Coord::unit() const
{
	Coord c;
	double length = this->abs();
	if(ALMOST_EQUAL_SMALL(length,0.0))
	{
//		fprintf(stdout,"Coord Error: unit, divided by zero\n");
		return *this;
	}
	c.pos[0] = pos[0]/length;
	c.pos[1] = pos[1]/length;
	c.pos[2] = pos[2]/length;
	return c;
}

inline bool Coord::normalize()
{
	double length = this->abs();
	if(ALMOST_EQUAL_SMALL(length,0.0))
	{
//		fprintf(stdout,"Coord Error: normalize, divided by zero\n");
		return false;
	}
	pos[0] /= length;
	pos[1] /= length;
	pos[2] /= length;
	return true;
}

// Transformer

// Change Coordinate as pos[2]->pos[0], pos[0]->pos[1], pos[1]->pos[2]
inline Coord Coord::Rotate(double delta_theta, double delta_phi, double theta)
{
	Coord v;
	v.pos[0] = this->pos[2];
	v.pos[1] = this->pos[0];
	v.pos[2] = this->pos[1];
	Coord nv;
	if(ALMOST_EQUAL_SMALL(delta_phi, 0.0))
	{		
		nv.pos[0] = v.pos[0]*cos(delta_theta)+v.pos[1]*sin(delta_theta);
		nv.pos[1] = v.pos[0]*(-1.0)*sin(delta_theta)+v.pos[1]*cos(delta_theta);
		nv.pos[2] = v.pos[2];
	}
	else
	{
		Coord Y(cos(theta),sin(theta),0.0);
		Coord Z(0.0,0.0,1.0);
		Coord X = cross(Y,Z);
		X.normalize();
		Coord pos;
		pos.pos[0] = dot(v,X);
		pos.pos[1] = dot(v,Y);
		pos.pos[2] = dot(v,Z);
		Coord newpos;
		newpos.pos[2] = pos.pos[2]*cos(delta_phi)+pos.pos[1]*sin(delta_phi);
		newpos.pos[1] = pos.pos[2]*(-1.0)*sin(delta_phi)+pos.pos[1]*cos(delta_phi);
		newpos.pos[0] = pos.pos[0];
		nv = X*newpos.pos[0]+Y*newpos.pos[1]+Z*newpos.pos[2];
	}

	return Coord(nv.pos[1],nv.pos[2],nv.pos[0]);
}

inline Coord Coord::Scale(Coord scale)
{
	Coord v = *this;
	for(int i = 0; i < 3; i ++)
		v[i] *= scale[i];
	return v;
}

inline Coord Coord::Spin(Coord SpinAxes, Coord Center, double SpinAngle)
{
	int i;
	Coord v = *this;
	if(SpinAxes.abs() < 1.0e-12)
		return v;
	
	SpinAxes.normalize();
	double d = -1.0*dot(SpinAxes,Center);
	for(i = 0; i < 3; i ++)
		if(!ALMOST_EQUAL_SMALL(SpinAxes[i], 0.0))
			break;
	Coord U(0.0,0.0,0.0);
	U[(i+2)%3] = Center[(i+2)%3]+1.0;
	U[(i+1)%3] = Center[(i+1)%3]+1.0;
	U[i] = (-1.0*d-SpinAxes[(i+2)%3]*U[(i+2)%3]-SpinAxes[(i+1)%3]*U[(i+1)%3])/SpinAxes[i];
	U -= Center;
	U.normalize();
	Coord V = cross(SpinAxes,U);
	V.normalize();
	Coord pos;
	pos.pos[0] = dot(v,U);
	pos.pos[1] = dot(v,V);
	pos.pos[2] = dot(v,SpinAxes);
	Coord newpos;
	newpos.pos[0] = pos.pos[0]*cos(SpinAngle)+pos.pos[1]*sin(SpinAngle);
	newpos.pos[1] = pos.pos[0]*(-1.0)*sin(SpinAngle)+pos.pos[1]*cos(SpinAngle);
	newpos.pos[2] = pos.pos[2];
	
	return U*newpos[0]+V*newpos[1]+SpinAxes*newpos[2];
}

inline double Coord::x() const
{
	return pos[0];
}
inline double Coord::y() const
{
	return pos[1];
}
inline double Coord::z() const
{
	return pos[2];
}

/* ================== Coordinate (2D) ================== */

// Initializer
inline void Coord2D::setCoords(double _x /* = 0.0 */, double _y /* = 0.0 */)
{
    pos[0] = _x;
    pos[1] = _y;
}

// Basic operator
inline Coord2D Coord2D::operator +(const Coord2D &c) const
{
	Coord2D coord;
	coord.pos[0] = pos[0]+c.pos[0];
	coord.pos[1] = pos[1]+c.pos[1];
	return coord;
}

inline Coord2D Coord2D::operator -(const Coord2D &c) const
{
	Coord2D coord;
	coord.pos[0] = pos[0]-c.pos[0];
	coord.pos[1] = pos[1]-c.pos[1];
	return coord;
}

inline Coord2D Coord2D::operator *(double v) const
{
	Coord2D coord;
	coord.pos[0] = pos[0]*v;
	coord.pos[1] = pos[1]*v;
	return coord;
}

inline Coord2D Coord2D::operator /(double v) const
{
	Coord2D coord;
	if(ALMOST_EQUAL_SMALL(v,0.0))
	{
		fprintf(stdout,"Coord2D Error: divided by zero\n");
		return coord;
	}
	coord.pos[0] = pos[0]/v;
	coord.pos[1] = pos[1]/v;
	return coord;
}

inline Coord2D Coord2D::operator -() const
{
	Coord2D coord;
	coord.pos[0] = -1.0*pos[0];
	coord.pos[1] = -1.0*pos[1];
	return coord;
}

//
inline Coord2D& Coord2D::operator +=(const Coord2D &c)
{
	pos[0] += c.pos[0];
	pos[1] += c.pos[1];
	return *this;
}

inline Coord2D& Coord2D::operator -=(const Coord2D &c)
{
	pos[0] -= c.pos[0];
	pos[1] -= c.pos[1];
	return *this;
}

inline Coord2D& Coord2D::operator *=(double v)
{
	pos[0] *= v;
	pos[1] *= v;
	return *this;
}

inline Coord2D& Coord2D::operator /=(double v)
{
	if(ALMOST_EQUAL_SMALL(v,0.0))
	{
		fprintf(stdout,"Coord2D Error: divided by zero\n");
		return *this;
	}
	pos[0] /= v;
	pos[1] /= v;
	return *this;
}

//
inline Coord2D& Coord2D::operator =(const Coord2D &c)
{
    pos[0] = c.pos[0];
    pos[1] = c.pos[1];
    return *this;
}

//
inline bool Coord2D::operator ==(const Coord2D &c) const
{
	if(pos[0]!=c.pos[0] || pos[1]!=c.pos[1])
		return false;
	return true;
}

inline bool Coord2D::operator !=(const Coord2D &c) const
{
	if(pos[0]!=c.pos[0] || pos[1]!=c.pos[1])
		return true;
	return false;
}

//
inline double& Coord2D::operator [](int i)
{
    return pos[i];
}

inline double Coord2D::operator [](int i) const
{
    return pos[i];
}

// Advanced operator
inline double dot(const Coord2D& c1, const Coord2D& c2)
{
	return (c1.pos[0]*c2.pos[0] + c1.pos[1]*c2.pos[1]);
}

inline double angle(const Coord2D& c1, const Coord2D& c2)
{
	Coord2D tmp_c1(c1), tmp_c2(c2);
	tmp_c1.normalize();
	tmp_c2.normalize();

	double dot_ = dot(tmp_c1, tmp_c2);
	if(dot_ < -1.0)
		dot_ = -1.0;
	else if(dot_ > 1.0)
		dot_ = 1.0;
	
	return acos(dot_);
}

//
inline double Coord2D::sqrabs() const
{
	return (pos[0]*pos[0] + pos[1]*pos[1]);
}

inline double Coord2D::abs() const
{
	return sqrt(pos[0]*pos[0] + pos[1]*pos[1]);
}

//
inline Coord2D Coord2D::unit() const
{
	Coord2D c;
	double length = this->abs();
	if(ALMOST_EQUAL_SMALL(length, 0.0))
	{
		fprintf(stdout,"Coord2D Error: unit, divided by zero\n");
		return c;
	}
	c.pos[0] = pos[0]/length;
	c.pos[1] = pos[1]/length;
	return c;
}

inline bool Coord2D::normalize()
{
	double length = this->abs();
	if(ALMOST_EQUAL_SMALL(length, 0.0))
	{
		fprintf(stdout,"Coord2D Error: normalize, divided by zero\n");
		return false;
	}
	pos[0] /= length;
	pos[1] /= length;
	return true;
}
inline double Coord2D::x() const
{
	return pos[0];
}
inline double Coord2D::y() const
{
	return pos[1];
}


/* ================== Color (RGB) ================== */
#include <stdlib.h>
#include <time.h>

// Initializer
inline void Color::setColor(int r, int g, int b)
{
	rgb[0] = r/255.0;
    rgb[1] = g/255.0;
    rgb[2] = b/255.0;

    Clamp();
}

inline void Color::setColor(double r, double g, double b)
{
	rgb[0] = r;
	rgb[1] = g;
	rgb[2] = b;

    Clamp();
}

inline void Color::setRandomColor()
{
	for(int i = 0; i < 3; i ++)
		rgb[i] = (double)(rand()%255)/(double)255;
//	fprintf(stdout,"Random Color = (%f,%f,%f)\n",rgb[0],rgb[1],rgb[2]);
}

// Basic operator
inline Color Color::operator +(const Color& color) const
{
	Color c;
	for(int i = 0; i < 3; i ++)
		c[i] = rgb[i]+color[i];
	return c;
}

inline Color Color::operator -(const Color& color) const
{
	Color c;
	for(int i = 0; i < 3; i ++)
		c[i] = rgb[i]-color[i];
	return c;
}

inline Color Color::operator *(double v) const
{
	Color c;
	for(int i = 0; i < 3; i ++)
		c[i] = rgb[i]*v;
	return c;
}

inline Color Color::operator /(double v) const
{
	Color c;
	for(int i = 0; i < 3; i ++)
		c[i] = rgb[i]/v;
	return c;
}

inline Color Color::operator -() const
{
	Color c;
	for(int i = 0; i < 3; i ++)
		c[i] = -1.0*rgb[i];
	return c;
}

//
inline Color& Color::operator +=(const Color &color)
{
	for(int i = 0; i < 3; i ++)
		rgb[i] += color[i];
	return *this;
}

inline Color& Color::operator -=(const Color &color)
{
	for(int i = 0; i < 3; i ++)
		rgb[i] -= color[i];
	return *this;
}

inline Color& Color::operator *=(double v)
{
	for(int i = 0; i < 3; i ++)
		rgb[i] *= v;
	return *this;
}

inline Color& Color::operator /=(double v)
{
	for(int i = 0; i < 3; i ++)
		rgb[i] /= v;
	return *this;
}

//
inline Color& Color::operator=(const Color &color)
{
	for(int i = 0; i < 3; i ++)
		rgb[i] = color[i];
	return *this;
}

//
inline bool Color::operator ==(const Color &color) const
{
	for(int i = 0; i < 3; i ++)
		if(rgb[i] != color[i])
			return false;
	return true;
}

inline bool Color::operator !=(const Color &color) const
{
	for(int i = 0; i < 3; i ++)
		if(rgb[i] != color[i])
			return true;
	return false;
}

inline void Color::Clamp()
{
    for(int i = 0; i < 3; ++ i)
    {
        if(rgb[i] < 0.0)
            rgb[i] = 0.0;
        else if(rgb[i] > 1.0)
            rgb[i] = 1.0;
    }
}

//inline COLORREF Color::Color2Ref()
//{
//	return RGB((BYTE)(rgb[0]*255), (BYTE)(rgb[1]*255), (BYTE)(rgb[2]*255));
//}
//
//Color::Color(COLORREF color)
//{
//	rgb[0] = (double)GetRValue(color)/(double)255;
//	rgb[1] = (double)GetGValue(color)/(double)255;
//	rgb[2] = (double)GetBValue(color)/(double)255;
//}


/* ================== Color (RGB) ================== */
inline CCurvature& CCurvature::operator = (const CCurvature &curv)
{
	this->m_direction_kmax = curv.m_direction_kmax;
	this->m_direction_kmin = curv.m_direction_kmin;
	this->m_kmax = curv.m_kmax;
	this->m_kmin = curv.m_kmin;

	return *this;
}
