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
// MeshModel.cpp
//
// [Developer]
// Xu Dong
// State Key Lab of CAD&CG, Zhejiang University
// 
// [Date]
// 2005-08-09
//
// [Goal]
// Defining the general interface for mesh model
#include "MeshModel.h"
using namespace std;
// Constructor
MeshModel::MeshModel()
{
    m_IO.AttachKernel(&m_Kernel);

    m_Render.AttachKernel(&m_Kernel);
    m_Render.AttachAuxData(&m_AuxData);

    m_BasicOp.AttachKernel(&m_Kernel);
    m_BasicOp.AttachAuxData(&m_AuxData);

    m_AdvancedOp.AttachKernel(&m_Kernel);
    m_AdvancedOp.AttachAuxData(&m_AuxData);
    m_AdvancedOp.AttachBasicOp(&m_BasicOp);

    ClearData();
}

// Destructor
MeshModel::~MeshModel()
{
    ClearData();
}

// Initialization
void MeshModel::ClearData()
{
    m_Kernel.ClearData();
    m_AuxData.ClearData();
    m_IO.ClearData();
    m_Render.ClearData();
    m_BasicOp.ClearData();
    m_AdvancedOp.ClearData();

    m_bAttachModel = false;
}

// Input/Output functions
void MeshModel::AttachModel(string filename)
{
    if(!m_IO.LoadModel(filename))
        return;

    m_bAttachModel = true;
    m_BasicOp.InitModel();
	m_ModelName = filename;
}
string MeshModel::GetModelFileName()
{
	return m_ModelName;
}
void MeshModel::StoreModel(string filename)
{
//	m_AdvancedOp.CompactModel();
    m_IO.StoreModel(filename);
}

// Create a model by coord array and face index array
void MeshModel::CreateModel(CoordArray CoordArr, PolyIndexArray FaceArr)
{
    size_t nVertex = CoordArr.size();
    size_t nFace = FaceArr.size();
    m_Kernel.GetVertexInfo().GetCoord() = CoordArr;
    m_Kernel.GetFaceInfo().GetIndex()   = FaceArr;

    m_bAttachModel = true;
    m_BasicOp.InitModel();
    
    printf("#Vertex = %d, #Face = %d\n\n", nVertex, nFace);
}

// Rendering
void MeshModel::DrawModel()
{
    m_Render.DrawModel();
}

void MeshModel::UpdateModel()
{
    m_BasicOp.CalFaceNormal();
    m_BasicOp.CalVertexNormal();
    m_BasicOp.CalBoundingBox();
}

// Compact model means that rearrange the vertex index and face index, 
// so that invalid vertices/faces are excluded.
void MeshModel::CompactModel()
{
    m_AdvancedOp.CompactModel();
    m_BasicOp.InitModel();
}

// Transformation functions
void MeshModel::Translate(Coord delta)
{
    m_AdvancedOp.Translate(delta);
}

void MeshModel::Rotate(Coord center, Coord axis, double angle)
{
    m_AdvancedOp.Rotate(center, axis, angle);
}

void MeshModel::Scale(Coord center, Coord ratio)
{
    m_AdvancedOp.Scale(center, ratio);
}

int MeshModel::CreateTexture(const std::string& file_name)
{
	return m_Render.CreateTexture(file_name);
}
