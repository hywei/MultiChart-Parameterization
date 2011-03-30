// 2005-08-05
//
// [Goal]
// A general, flexible, versatile and easy-to-use mesh library for research purpose.
// Supporting arbitrary polygonal meshes as input, but with a 
// primary focus on triangle meshes.

/* ================== File Information ================== */
// [Name]
// MeshModelAuxData.h
//
// [Developer]
// Xu Dong
// State Key Lab of CAD&CG, Zhejiang University
// 
// [Date]
// 2005-08-09
//
// [Goal]
// Defining the auxiliary database for mesh model
// For the purpose of debugging and visualization



#include "../Common/Utility.h"
#include <vector>
#include <string>
#pragma once



/* ================== Auxiliary Element - Point Information ================== */

class PointInfo
{

private:
    // The following properties (arrays) are one-for-each-point
    CoordArray      m_Coord;    // Point coordinate array
    ColorArray      m_Color;    // Point color array
    Utility         util;

public:
    // Constructor
    PointInfo();

    // Destructor
    ~PointInfo();

    // Initializer
    void ClearData();

    // Get/Set functions
    CoordArray& GetCoord() { return m_Coord; }
    ColorArray& GetColor() { return m_Color; }

    // Add/Remove functions
    void Add(Coord p, Color c);
    void Remove(int idx);

};



/* ================== Auxiliary Element - Line Information ================== */

class LineInfo
{
private:
    // The following properties (arrays) are one-for-each-line
    CoordArray      m_fCoord;   // First point coordinate array
    CoordArray      m_sCoord;   // Second point coordinate array
    ColorArray      m_Color;    // Line color array
    
public:
    // Constructor
    LineInfo();

    // Destructor
    ~LineInfo();

    // Initializer
    void ClearData();

    // Get/Set functions
    CoordArray& GetfCoord() { return m_fCoord; }
    CoordArray& GetsCoord() { return m_sCoord; }
    ColorArray& GetColor() { return m_Color; }

    // Add/Remove functions
    void Add(Coord fp, Coord sp, Color c);
    void Remove(int idx);

};



/* ================== Auxiliary Element - Polygon Information ================== */

class PolygonInfo
{
private:
    // The following properties (arrays) are one-for-each-polygon
    PolyCoordArray  m_Coords;   // Polygon coordinate array
    ColorArray      m_Color;    // Polygon color array

public:
    // Constructor
    PolygonInfo();

    // Destructor
    ~PolygonInfo();

    // Initializer
    void ClearData();

    // Get/Set functions
    PolyCoordArray& GetCoords() { return m_Coords; }
    ColorArray& GetColor() { return m_Color; }

    // Add/Remove functions
    void Add(CoordArray& vCoord, Color c);
	void Add(CoordArray& vCoord, ColorArray& cArray);
    void Remove(int idx);

};



/* ================== Auxiliary Data Information ================== */

class MeshModelAuxData
{
private:
    // Basic Mesh Model Kernel = Vertex Coord Array + Face Index Array
    PointInfo   m_PointInfo;
    LineInfo    m_LineInfo;
    PolygonInfo m_PolygonInfo;

public:
	// Constructor & Destructor
	MeshModelAuxData();
	virtual ~MeshModelAuxData();

    // Initializer
    void ClearData();

    // Get functions
    PointInfo&      GetPointInfo()   { return m_PointInfo; }
    LineInfo&       GetLineInfo()    { return m_LineInfo; }
    PolygonInfo&    GetPolygonInfo() { return m_PolygonInfo; }

    // Add/Remove functions
    void AddPoint(Coord p, Color c);
    void AddLine(Coord fp, Coord sp, Color c);
    void AddPolygon(CoordArray& vCoord, Color c);

    void RemovePoint(int idx);
    void RemoveLine(int idx);
    void RemovePolygon(int idx);

};
