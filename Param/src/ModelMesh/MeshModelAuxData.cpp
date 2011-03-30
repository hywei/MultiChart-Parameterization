// 2005-08-05
//
// [Goal]
// A general, flexible, versatile and easy-to-use mesh library for research purpose.
// Supporting arbitrary polygonal meshes as input, but with a 
// primary focus on triangle meshes.

/* ================== File Information ================== */
// [Name]
// MeshModelAuxData.cpp
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

#include "MeshModelAuxData.h"



/* ================== Auxiliary Element - Point Information ================== */

// Constructor
PointInfo::PointInfo()
{

}

// Destructor
PointInfo::~PointInfo()
{

}

// Initializer
void PointInfo::ClearData()
{
    util.FreeVector(m_Coord);
    util.FreeVector(m_Color);
}

// Add/Remove functions
void PointInfo::Add(Coord p, Color c)
{
    m_Coord.push_back(p);
    m_Color.push_back(c);
}

void PointInfo::Remove(int idx)
{
    int n = (int) m_Coord.size();
    if(idx >= 0 && idx < n)
    {
        m_Coord.erase(m_Coord.begin()+idx);
        m_Color.erase(m_Color.begin()+idx);
    }
}



/* ================== Auxiliary Element - Line Information ================== */

// Constructor
LineInfo::LineInfo()
{

}

// Destructor
LineInfo::~LineInfo()
{

}

// Initializer
void LineInfo::ClearData()
{
    Utility util;
    util.FreeVector(m_fCoord);
    util.FreeVector(m_sCoord);
    util.FreeVector(m_Color);
}

// Add/Remove functions
void LineInfo::Add(Coord fp, Coord sp, Color c)
{
    m_fCoord.push_back(fp);
    m_sCoord.push_back(sp);
    m_Color.push_back(c);
}

void LineInfo::Remove(int idx)
{
    int n = (int) m_fCoord.size();
    if(idx >= 0 && idx < n)
    {
        m_fCoord.erase(m_fCoord.begin()+idx);
        m_sCoord.erase(m_sCoord.begin()+idx);
        m_Color.erase(m_Color.begin()+idx);
    }
}



/* ================== Auxiliary Element - Polygon Information ================== */

// Constructor
PolygonInfo::PolygonInfo()
{

}

// Destructor
PolygonInfo::~PolygonInfo()
{

}

// Initializer
void PolygonInfo::ClearData()
{
    Utility util;
    util.FreeVector(m_Coords);
    util.FreeVector(m_Color);
}

// Add/Remove functions
void PolygonInfo::Add(CoordArray& vCoord, Color c)
{
    m_Coords.push_back(vCoord);
    m_Color.push_back(c);
}

void PolygonInfo::Remove(int idx)
{
    int n = (int) m_Coords.size();
    if(idx >= 0 && idx < n)
    {
        m_Coords.erase(m_Coords.begin()+idx);
        m_Color.erase(m_Color.begin()+idx);
    }
}



/* ================== Auxiliary Data Information ================== */



// Constructor
MeshModelAuxData::MeshModelAuxData()
{
    
}

// Destructor
MeshModelAuxData::~MeshModelAuxData()
{

}

// Initializer
void MeshModelAuxData::ClearData()
{
    m_PointInfo.ClearData();
    m_LineInfo.ClearData();
    m_PolygonInfo.ClearData();
}

// Add/Remove functions
void MeshModelAuxData::AddPoint(Coord p, Color c)
{
    m_PointInfo.Add(p, c);
}

void MeshModelAuxData::AddLine(Coord fp, Coord sp, Color c)
{
    m_LineInfo.Add(fp, sp, c);
}

void MeshModelAuxData::AddPolygon(CoordArray& vCoord, Color c)
{
    m_PolygonInfo.Add(vCoord, c);
}

void MeshModelAuxData::RemovePoint(int idx)
{
    m_PointInfo.Remove(idx);
}

void MeshModelAuxData::RemoveLine(int idx)
{
    m_LineInfo.Remove(idx);
}

void MeshModelAuxData::RemovePolygon(int idx)
{
    m_PolygonInfo.Remove(idx);
}