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
// MeshModel.h
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


#include "BasicElement.h"       // Basic elements -- including Edge, Face, Point, Line, etc

#include "../Common/Utility.h"    // Utility support

#include "MeshModelKernel.h"    // Kernel database and states
#include "MeshModelAuxData.h"   // Auxiliary database and states -- for debugging and visualization
#include "MeshModelIO.h"        // Input/output functions
#include "MeshModelRender.h"    // Rendering functions
#include "MeshModelBasicOp.h"       // Basic operations -- those do not change mesh topology
#include "MeshModelAdvancedOp.h"    // Advanced operations -- those change mesh topology

#pragma once



/* ================== General Mesh Interface ================== */

class MeshModel
{
public:
    MeshModelKernel     m_Kernel;
    MeshModelAuxData    m_AuxData;
    MeshModelIO         m_IO;
    MeshModelRender     m_Render;
    MeshModelBasicOp    m_BasicOp;
    MeshModelAdvancedOp m_AdvancedOp;
    bool        m_bAttachModel;
	std::string      m_ModelName;

private:
    Utility     m_Utility;

public:
    // Constructor
    MeshModel();

    // Destructor
    ~MeshModel();

    // Initialization
    void ClearData();

    // Input/Output functions
    void AttachModel(std::string filename);
    void StoreModel(std::string filename);

    // Create a model by coord array and face index array
    void CreateModel(CoordArray CoordArr, PolyIndexArray FaceArr);
    std::string GetModelFileName();

    // Rendering
    void DrawModel();

    // Update model
    void UpdateModel();

    // Compact model
    void CompactModel();

    // Transformation functions
    void Translate(Coord delta);
    void Rotate(Coord center, Coord axis, double angle);
    void Scale(Coord center, Coord ratio);
    

	int CreateTexture(const std::string& file_name);
};
