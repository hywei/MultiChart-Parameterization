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
// MeshModelKernel.cpp
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

#include "MeshModelKernel.h"



/* ================== Kernel Element - Mesh Vertex Information ================== */

// Constructor
VertexInfo::VertexInfo()
{

}

// Destructor
VertexInfo::~VertexInfo()
{

}

// Initializer
void VertexInfo::ClearData()
{
    Utility util;
    util.FreeVector(m_Coord);
    util.FreeVector(m_Normal);
    util.FreeVector(m_Color);
    util.FreeVector(m_TexCoord);
    util.FreeVector(m_Flag);

    util.FreeVector(m_AdjFaces);
    util.FreeVector(m_AdjVertices);
    util.FreeVector(m_AdjEdges);

    m_nVertices = 0;
}



/* ================== Kernel Element - Mesh Face Information ================== */

// Constructor
FaceInfo::FaceInfo()
{

}

// Destructor
FaceInfo::~FaceInfo()
{

}

// Initializer
void FaceInfo::ClearData()
{
    Utility util;
    util.FreeVector(m_Index);
    util.FreeVector(m_Normal);
    util.FreeVector(m_Color);
    util.FreeVector(m_TexCoord);
    util.FreeVector(m_Flag);
	util.FreeVector(m_FaceBaryCenter);
	util.FreeVector(m_TexIndex);

    m_nFaces = 0;
}



/* ================== Kernel Element - Mesh Edge Information ================== */

// Constructor
EdgeInfo::EdgeInfo()
{

}

// Destructor
EdgeInfo::~EdgeInfo()
{

}

// Initializer
void EdgeInfo::ClearData()
{
    Utility util;
    util.FreeVector(m_VtxIndex);
    util.FreeVector(m_Color);
    util.FreeVector(m_Flag);
	util.FreeVector(m_DihedralAngle);
	util.FreeVector(m_FaceIndex);

    m_nHalfEdges = 0;
}



/* ================== Kernel Element - Mesh Model Information ================== */

// Constructor
ModelInfo::ModelInfo()
{

}

// Destructor
ModelInfo::~ModelInfo()
{

}

// Initializer
void ModelInfo::ClearData()
{
    m_Flag = 0;
    m_FileName.resize(0);

    m_Type = MODEL_TYPE_POLYGON_SOAP;
    
    m_nVertices = m_nFaces = m_nHalfEdges = 0;

    m_nBoundaries = 0;
    m_nComponents = 0;

    m_BoxMin = m_BoxMax = m_BoxDim = m_SphereCenter = Coord(0.0, 0.0, 0.0);
    m_SphereRadius = 0.0;

    m_AvgEdgeLength = 0.0;
    m_AvgFaceArea = 0.0;

    util.FreeVector(m_Boundaries);
}

void ModelInfo::SetFileName(std::string filename)
{
    m_FileName = filename;
}

bool ModelInfo::IsTriMesh()
{
    return util.IsSetFlag(m_Flag, MODEL_FLAG_TRIMESH);
}

bool ModelInfo::IsQuadMesh()
{
    return util.IsSetFlag(m_Flag, MODEL_FLAG_QUADMESH);
}

bool ModelInfo::IsGeneralMesh()
{
    return util.IsSetFlag(m_Flag, MODEL_FLAG_GENERALMESH);
}

bool ModelInfo::IsClosed()
{
    return (util.IsSetFlag(m_Flag, MODEL_FLAG_MANIFOLD) && m_nBoundaries == 0);
}

bool ModelInfo::IsManifold()
{
    return util.IsSetFlag(m_Flag, MODEL_FLAG_MANIFOLD);
}

bool ModelInfo::IsTriManifold()
{
    return (util.IsSetFlag(m_Flag, MODEL_FLAG_MANIFOLD) && util.IsSetFlag(m_Flag, MODEL_FLAG_TRIMESH));
}

bool ModelInfo::IsPatch()
{
    return (util.IsSetFlag(m_Flag, MODEL_FLAG_MANIFOLD) && m_nBoundaries == 1);
}

/* ================== Mesh Model Kernel  ================== */

// Constructor
MeshModelKernel::MeshModelKernel()
{
    
}

// Destructor
MeshModelKernel::~MeshModelKernel()
{

}

// Initializer
void MeshModelKernel::ClearData()
{
    m_VertexInfo.ClearData();
    m_FaceInfo.ClearData();
    m_EdgeInfo.ClearData();
    m_ModelInfo.ClearData();
}
