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
// MeshModelAdvancedOp.h
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



#include "MeshModelKernel.h"
#include "MeshModelAuxData.h"
#include "MeshModelBasicOp.h"
#include "../Common/Utility.h"

#pragma once



/* ================== Advanced Operation ================== */

class MeshModelAdvancedOp
{
private:
    MeshModelKernel* kernel;
    MeshModelAuxData* auxdata;
    MeshModelBasicOp* basicop;
    Utility util;

public:
    // Constructor
    MeshModelAdvancedOp();

    // Destructor
    ~MeshModelAdvancedOp();

    // Initializer
    void ClearData();
    void AttachKernel(MeshModelKernel* pKernel);
    void AttachAuxData(MeshModelAuxData* pAuxData);
    void AttachBasicOp(MeshModelBasicOp* pBasicOp);

    // Compact model
    void CompactModel();

    // Euler operator
    void EdgeCollapse(VertexID vID1, VertexID vID2);
    void VertexSplit(VertexID vID);
    void EdgeSplit(VertexID vID1, VertexID vID2);

    // Model transformation
    void Translate(Coord delta);
    void Rotate(Coord center, Coord axis, double angle);
    void Scale(Coord center, Coord ratio);

};