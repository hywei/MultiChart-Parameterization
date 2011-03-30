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
// MeshModelKernel.h
//
// [Developer]
// Xu Dong
// State Key Lab of CAD&CG, Zhejiang University
// 
// [Date]
// 2005-08-05
//
// [Goal]
// Defining the kernel components of the mesh model
//
// This class includes four kernel components - vertex, face, edge(optional) and model components.
// This class is responsible for the kernel data storage and component state information



#include "../Common/Utility.h"
#include <vector>
#include <string>
#include <iostream>

#pragma once


/* ================== Kernel Element Flags ================== */
// General Flags
#define FLAG_INVALID        0X00000001
#define FLAG_HIGHLIGHTED    0X00000002
#define FLAG_SELECTED       0X00000004
#define FLAG_VISITED        0X00000008

#define FLAG_USER1          0X00000010
#define FLAG_USER2          0X00000020
#define FLAG_USER3          0X00000040
#define FLAG_USER4          0X00000080

// Vertex Flags
#define VERTEX_FLAG_ISOLATED    0X00000100
#define VERTEX_FLAG_MANIFOLD    0X00000200
#define VERTEX_FLAG_BOUNDARY    0X00000400

// Face Flags
#define FACE_FLAG_BOUNDARY      0X00000100
#define FACE_FLAG_MANIFOLD      0X00000200

// Edge Flags
#define EDGE_FLAG_BOUNDARY      0X00000100
#define EDGE_FLAG_MANIFOLD      0X00000200

// Model Flags
#define MODEL_FLAG_TRIMESH      0X00000100
#define MODEL_FLAG_QUADMESH     0X00000200
#define MODEL_FLAG_GENERALMESH  0X00000400

#define MODEL_FLAG_MANIFOLD     0X00000800

// Model Type
#define MODEL_TYPE_POINT_CLOUD  0X00000001
#define MODEL_TYPE_POLYGON_SOAP 0X00000002



/* ================== Kernel Element - Mesh Vertex Information ================== */

class VertexInfo
{
private:
    // The following properties (arrays) are one-for-each-vertex
    CoordArray      m_Coord;    // Vertex coordinate array
    NormalArray     m_Normal;   // Vertex normal array
    ColorArray      m_Color;    // Vertex color array
    TexCoordArray   m_TexCoord; // Vertex texture coordinate array
    FlagArray       m_Flag;     // Vertex 32-bit flag array

    PolyIndexArray  m_AdjFaces;     // Vertex adjacent face-index array
    PolyIndexArray  m_AdjVertices;  // Vertex adjacent vertex-index array
    PolyIndexArray  m_AdjEdges;     // Vertex adjacent half-edge-index array   
    
	CurvatureArray  m_Curvatures;   // Vertex curvature array

    int m_nVertices;

public:
    // Constructor
    VertexInfo();

    // Destructor
    ~VertexInfo();

    // Initializer
    void ClearData();

    // Get/Set functions
    CoordArray& GetCoord() { return m_Coord; }
    NormalArray& GetNormal() { return m_Normal; }
    ColorArray& GetColor() { return m_Color; }
    TexCoordArray& GetTexCoord() { return m_TexCoord; }
    FlagArray& GetFlag() { return m_Flag; }

    PolyIndexArray& GetAdjFaces() { return m_AdjFaces; }
    PolyIndexArray& GetAdjVertices() { return m_AdjVertices; }
    PolyIndexArray& GetAdjEdges() { return m_AdjEdges; }
	CurvatureArray& GetCurvatures() { return m_Curvatures; }
};



/* ================== Kernel Element - Mesh Face Information ================== */

class FaceInfo
{
private:
    // The following properties (arrays) are one-for-each-face
    PolyIndexArray     m_Index;    // Face vertex-index array
    NormalArray        m_Normal;   // Face normal array
    ColorArray         m_Color;    // Face color array
    PolyTexCoordArray  m_TexCoord; // Face vertex-texture-coordinate array
    FlagArray          m_Flag;     // Face 32-bit flag array
	CoordArray         m_FaceBaryCenter; // the barycenter of each face.
	DoubleArray        m_FaceArea;	
	PolyIndexArray     m_TexIndex; // Face vertex-texture-index array
     
    int m_nFaces;

public:
    // Constructor
    FaceInfo();

    // Destructor
    ~FaceInfo();

    // Initializer
    void ClearData();

    // Get/Set functions
    PolyIndexArray& GetIndex() { return m_Index; }
    NormalArray& GetNormal() { return m_Normal; }
    ColorArray& GetColor() { return m_Color; }
    PolyTexCoordArray& GetTexCoord() { return m_TexCoord; }
    FlagArray& GetFlag() { return m_Flag; }
	CoordArray& GetBaryCenter() { return m_FaceBaryCenter; }
	DoubleArray& GetFaceArea() { return m_FaceArea; }
	PolyIndexArray& GetTexIndex() { return m_TexIndex; }

};



/* ================== Kernel Element - Mesh Edge Information ================== */

class EdgeInfo
{
private:
    PolyIndexArray  m_VtxIndex;    // Edge vertex-index array
    ColorArray      m_Color;    // Edge color array
    FlagArray       m_Flag;     // Edge 32-bit flag array
    DoubleArray     m_DihedralAngle;     // the dihedral angle of each edge.
	PolyIndexArray  m_FaceIndex;    // Edge face-index array
    int m_nHalfEdges;
    
public:
    // Constructor
    EdgeInfo();

    // Destructor
    ~EdgeInfo();

    // Initializer
    void ClearData();

    // Get/Set functions
    PolyIndexArray& GetVertexIndex() { return m_VtxIndex; }
    ColorArray& GetColor() { return m_Color; }
    FlagArray& GetFlag() { return m_Flag; }
	DoubleArray& GetDihedralAngle() {return m_DihedralAngle;}
	PolyIndexArray& GetFaceIndex() {return m_FaceIndex;}

};



/* ================== Kernel Element - Mesh Model Information ================== */

class ModelInfo
{
private:
    Flag    m_Flag;     // Model 32-bit flag
	std::string  m_FileName; // File name of the mesh model (Path + Tiltle + Ext)
    int     m_Type;     // Model type - Point cloud or Polygon soap
    Color   m_Color;    // Model color

    int     m_nVertices;        // Number of valid vertices
    int     m_nFaces;           // Number of valid faces
    int     m_nHalfEdges;       // Number of valid half edges

    int     m_nBoundaries;  // Number of boundaries
    int     m_nComponents;  // Number of connected components

    Coord   m_BoxMin, m_BoxMax, m_BoxDim;
    Coord   m_SphereCenter;
    double  m_SphereRadius;

    double  m_AvgEdgeLength;
    double  m_AvgFaceArea;

    PolyIndexArray  m_Boundaries;   // Boundary vertex loop
    Utility     util;

public:
    // Constructor
    ModelInfo();

    // Destructor
    ~ModelInfo();

    // Initializer
    void ClearData();

    // Queries
    Flag& GetFlag() { return m_Flag; }

    std::string GetFileName() { return m_FileName; } 
    void SetFileName(std::string filename);
    
    int GetType() { return m_Type; }
    void SetType(int type) { m_Type = type; }

    Color GetColor() { return m_Color; }
    void SetColor(Color color) { m_Color = color; }
    
    int GetVertexNum() { return m_nVertices; }
    void SetVertexNum(int num) { m_nVertices = num; }

    int GetFaceNum() { return m_nFaces; }
    void SetFaceNum(int num) { m_nFaces = num; }

    int GetHalfEdgeNum() { return m_nHalfEdges; }
    void SetHalfEdgeNum(int num) { m_nHalfEdges = num; }
    
    int GetBoundaryNum() { return m_nBoundaries; }
    void SetBoundaryNum(int num) { m_nBoundaries = num; }

    int GetComponentNum() { return m_nComponents; }
    void SetComponentNum(int num) { m_nComponents = num; }

    void GetBoundingBox(Coord& BoxMin, Coord& BoxMax, Coord& BoxDim) { BoxMin = m_BoxMin; BoxMax = m_BoxMax; BoxDim = m_BoxDim; }
    void SetBoundingBox(Coord BoxMin, Coord BoxMax, Coord BoxDim) { m_BoxMin = BoxMin; m_BoxMax = BoxMax; m_BoxDim = BoxDim; }

    void GetBoundingSphere(Coord& Center, double& radius) { Center = m_SphereCenter; radius = m_SphereRadius; }
    void SetBoundingSphere(Coord Center, double Radius) { m_SphereCenter = Center; m_SphereRadius = Radius; }

    double GetAvgEdgeLength() { return m_AvgEdgeLength; }
    void SetAvgEdgeLength(double length) { m_AvgEdgeLength = length; }

    double GetAvgFaceArea() { return m_AvgFaceArea; }
    void SetAvgFaceArea(double area) { m_AvgFaceArea = area; }    

    PolyIndexArray& GetBoundary() { return m_Boundaries; }
    
    // Predictions    
    bool IsTriMesh();       // Whether the model is a triangle mesh (only containing triangles)
    bool IsQuadMesh();      // Whether the model is a quadrangle mesh (only containing quadangles)
    bool IsGeneralMesh();   // Whether the model is a general mesh (the rest, not a tri- or quad- one)

    bool IsClosed();        // Whether the model is a closed one (no boundaries)

    bool IsManifold();      // Whether the model is a 2-manifold one (locally disc-like)
    bool IsTriManifold();   // Whether the model is a 2-mainfold triangle mesh
    bool IsPatch();         // Whether the model is a patch (single boundary, 2-manifold)
};



/* ================== Mesh Model Kernel  ================== */

class MeshModelKernel
{
private:
    // Basic Mesh Model Kernel = Vertex Coord Array + Face Index Array
    VertexInfo  m_VertexInfo;
    FaceInfo    m_FaceInfo;
    EdgeInfo    m_EdgeInfo;
    ModelInfo   m_ModelInfo;
    Utility     util;

public:
	// Constructor & Destructor
	MeshModelKernel();
	virtual ~MeshModelKernel();

    // Initializer
    void ClearData();

    // Get functions
    VertexInfo& GetVertexInfo() { return m_VertexInfo; }
    EdgeInfo&   GetEdgeInfo()   { return m_EdgeInfo; }
    FaceInfo&   GetFaceInfo()   { return m_FaceInfo; }
    ModelInfo&  GetModelInfo()  { return m_ModelInfo; }
};
