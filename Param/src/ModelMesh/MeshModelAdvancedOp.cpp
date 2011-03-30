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
// MeshModelAdvancedOp.cpp
//
// [Developer]
// Xu Dong
// State Key Lab of CAD&CG, Zhejiang University
// 
// [Date]
// 2005-08-09
//
// [Goal]
// Defining the basic operations of the kernel components of the mesh model
// Including predictions, queries, Euler operations, etc
//
// This class is responsible for the manipulations of the kernel components and states
// Those change the topology of the mesh model


#include "MeshModelAdvancedOp.h"
#include "../Numerical/Rotation.h"
#include <cassert>


// Constructor
MeshModelAdvancedOp::MeshModelAdvancedOp()
{
    kernel = NULL ;
}

// Destructor
MeshModelAdvancedOp::~MeshModelAdvancedOp()
{

}

// Initializer
void MeshModelAdvancedOp::ClearData()
{
    
}

void MeshModelAdvancedOp::AttachKernel(MeshModelKernel* pKernel)
{
    assert(pKernel != NULL);
    kernel = pKernel;
}

void MeshModelAdvancedOp::AttachAuxData(MeshModelAuxData* pAuxData)
{
    assert(pAuxData != NULL);
    auxdata = pAuxData;
}

void MeshModelAdvancedOp::AttachBasicOp(MeshModelBasicOp* pBasicOp)
{
    assert(pBasicOp != NULL);
    basicop = pBasicOp;
}

// Compact model
void MeshModelAdvancedOp::CompactModel()
{
    CoordArray& vCoord = kernel->GetVertexInfo().GetCoord();
    FlagArray& vFlag = kernel->GetVertexInfo().GetFlag();
    PolyIndexArray& fIndex = kernel->GetFaceInfo().GetIndex();
    FlagArray& fFlag = kernel->GetFaceInfo().GetFlag();
    size_t nVertex = vCoord.size();
	size_t nFace = fIndex.size();

    IndexArray VtxIDMap, FaceIDMap;
    VtxIDMap.resize(nVertex);
    FaceIDMap.resize(nFace);

    size_t i;
    int index = 0;
    size_t nNewVtx, nNewFace;
    for(i = 0; i < nVertex; ++ i)
    {
        Flag vf = vFlag[i];
        if(util.IsSetFlag(vf, FLAG_INVALID))
        {
            VtxIDMap[i] = -1;
            printf("Invalid Vertex = %d\n", i); // Debug
        }
        else
        {
            VtxIDMap[i] = index ++;
        }
    }
    nNewVtx = index;

    index = 0;
    for(i = 0; i < nFace; ++ i)
    {
        Flag ff = fFlag[i];
        if(util.IsSetFlag(ff, FLAG_INVALID))
        {
            FaceIDMap[i] = -1;
        }
        else
        {
            FaceIDMap[i] = index ++;
        }
    }
    nNewFace = index;

    for(i = 0; i < nVertex; ++ i)
    {
        int idx = VtxIDMap[i];
        if(idx != -1 && idx != i)   // Invalid and avoid duplicate copy
            vCoord[idx] = vCoord[i];
    }
    vCoord.erase(vCoord.begin()+nNewVtx, vCoord.end());

    for(i = 0; i < nFace; ++ i)
    {
        int idx = FaceIDMap[i];
        if(idx == -1)
            continue;

        IndexArray& newf = fIndex[idx];
        IndexArray oldf = fIndex[i];
        size_t n = oldf.size();
        newf.resize(n);
        for(size_t j = 0; j < n; ++ j)
        {
            VertexID old_vID = oldf[j]; // Debug
            VertexID vID = VtxIDMap[oldf[j]];
            assert(vID != -1);
            newf[j] = vID;
        }
    }
    fIndex.erase(fIndex.begin()+nNewFace, fIndex.end());
}

// Model transformation
void MeshModelAdvancedOp::Translate(Coord delta)
{
    // Update vertex coords
    CoordArray& vCoord = kernel->GetVertexInfo().GetCoord();
    size_t i, n = vCoord.size();
    for(i = 0; i < n; ++ i)
        vCoord[i] += delta;

    // Update bounding box
    Coord BoxMin, BoxMax, BoxDim, Center;
    double radius;
    kernel->GetModelInfo().GetBoundingBox(BoxMin, BoxMax, BoxDim);
    kernel->GetModelInfo().GetBoundingSphere(Center, radius);

    BoxMin += delta;
    BoxMax += delta;
    Center += delta;
    kernel->GetModelInfo().SetBoundingBox(BoxMin, BoxMax, BoxDim);
    kernel->GetModelInfo().SetBoundingSphere(Center, radius);
}

// angle - degree
void MeshModelAdvancedOp::Rotate(Coord center, Coord axis, double angle)
{
    ROTATION rot(axis, angle);

    // Update vertex coords
    CoordArray& vCoord = kernel->GetVertexInfo().GetCoord();
    size_t i, n = vCoord.size();
    for(i = 0; i < n; ++ i)
        vCoord[i] = VectorTransform(vCoord[i]-center, rot) + center;

    // Update model
    basicop->CalFaceNormal();
    basicop->CalVertexNormal();
    basicop->CalBoundingBox();
}

void MeshModelAdvancedOp::Scale(Coord center, Coord ratio)
{
    // Update vertex coords
    CoordArray& vCoord = kernel->GetVertexInfo().GetCoord();
    int i, j, n = (int) vCoord.size();
    for(i = 0; i < n; ++ i)
        for(j = 0; j < 3; ++ j)
            vCoord[i][j] = (vCoord[i][j]-center[j])*ratio[j] + center[j];

    // Update bounding box
    Coord BoxMin, BoxMax, BoxDim, Center;
    double radius;
    kernel->GetModelInfo().GetBoundingBox(BoxMin, BoxMax, BoxDim);
    kernel->GetModelInfo().GetBoundingSphere(Center, radius);

    for(j = 0; j < 3; ++ j)
    {
        BoxMin[j] = (BoxMin[j]-center[j])*ratio[j] + center[j];
        BoxMax[j] = (BoxMax[j]-center[j])*ratio[j] + center[j];
        Center[j] = (Center[j]-center[j])*ratio[j] + center[j];
    }

    double r = std::max(ratio[0], ratio[1]);
    r = std::max(r, ratio[2]);
    radius *= r;
    kernel->GetModelInfo().SetBoundingBox(BoxMin, BoxMax, BoxDim);
    kernel->GetModelInfo().SetBoundingSphere(Center, radius);
}
