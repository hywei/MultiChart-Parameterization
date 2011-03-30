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
// MeshModelRender.h
//
// [Developer]
// Xu Dong
// State Key Lab of CAD&CG, Zhejiang University
// 
// [Date]
// 2005-08-05
//
// [Goal]
// Defining various rendering modes
// Including polygon (flag/smooth), wireframe, point, texture mapping, etc



#include "MeshModelKernel.h"
#include "MeshModelAuxData.h"
#include "../Common/Utility.h"
#include "../OpenGL/GLElement.h"
#pragma once



/* ================== Render Model Macros ================== */

#define RENDER_MODEL_VERTICES               0X00000001
#define RENDER_MODEL_WIREFRAME              0X00000002
#define RENDER_MODEL_SOLID_FLAT             0X00000003
#define RENDER_MODEL_SOLID_SMOOTH           0X00000004
#define RENDER_MODEL_SOLID_WITH_SHARPNESS   0X00000005
#define RENDER_MODEL_SOLID_AND_WIREFRAME    0X00000006
#define RENDER_MODEL_HIGHLIGHT_ONLY         0X00000007
#define RENDER_MODEL_TEXTURE_MAPPING        0X00000008
#define RENDER_MODEL_VERTEX_COLOR           0X00000009
#define RENDER_MODEL_FACE_COLOR             0X0000000A
#define RENDER_MODEL_TRANSPARENT            0X0000000B

#define RENDER_MODEL_AXIS                   0X00000100
#define RENDER_MODEL_BOUNDING_BOX           0X00000200
#define RENDER_MODEL_LIGHTING               0X00000400
#define RENDER_MODEL_COLOR                  0X00000800      //

#define RENDER_MODEL_POINT                  0X00010000
#define RENDER_MODEL_LINE                   0X00020000
#define RENDER_MODEL_POLYGON                0X00040000
#define RENDER_MODEL_ELEMENT                0X00080000

#define RENDER_MODEL_KMIN_CURVATURE         0X00100000
#define RENDER_MODEL_KMAX_CURVATURE         0X00200000

#define RENDER_MODEL_PARAM_TEXTURE          0X00400000

/* ================== Material ================== */

class Material
{
public:
	float m_Ambient[4];
	float m_Diffuse[4];
	float m_Specular[4];
	float m_Emissive[4];
	float m_Shininess;
	float m_Transparency;

public:
    // Constructor
    Material() {}
    
    Material(float ambient[4], float diffuse[4], float specular[4], float emissive[4], float shininess, float transparency = 0.0f)
    {
        for(int i = 0; i < 4; ++ i)
        {
            m_Ambient[i] = ambient[i];
            m_Diffuse[i] = diffuse[i];
            m_Specular[i] = specular[i];
            m_Emissive[i] = emissive[i];
            m_Shininess = shininess;
            m_Transparency = transparency;
        }
    }

    // Destructor
    ~Material() {}

    // Operator
    Material& operator =(Material& mat)
    {
        for(int i = 0; i < 4; ++ i)
        {
            m_Ambient[i] = mat.m_Ambient[i];
            m_Diffuse[i] = mat.m_Diffuse[i];
            m_Specular[i] = mat.m_Specular[i];
            m_Emissive[i] = mat.m_Emissive[i];
            m_Shininess = mat.m_Shininess;
            m_Transparency = mat.m_Transparency;
        }
        return *this;
    }
};



/* ================== Mesh Model Render ================== */

class MeshModelRender
{
private:
    int     m_Mode;
    int     m_State;
    
    Color   m_BbColor;
    
    Color   m_DftVtxColor;
    Color   m_DftEdgeColor;
    Color   m_DftFaceColor;
    
    Color   m_DftPointColor;
    Color   m_DftLineColor;
    Color   m_DftPolygonColor;
    Color   m_DftElementColor;

    Color   m_DftHltVtxColor;
    Color   m_DftHltEdgeColor;
    Color   m_DftHltFaceColor;

    float   m_DftVtxSize;
    float   m_DftEdgeWidth;
    
    float   m_DftPointSize;
    float   m_DftLineWidth;

    float   m_DftBdyVtxSize;
    float   m_DftBdyEdgeWidth;

    float   m_DftHltVtxSize;
    float   m_DftHltEdgeWidth;

    Material m_ModelMaterial;
    Material m_ElementMaterial;

	GLMaterial m_DefaultMaterial;

    MeshModelKernel* kernel;
    MeshModelAuxData* auxdata;

private:
    Utility util;

    // Various rendering functions
	void DrawModelVertices();
	void DrawModelWireframe();
	void DrawModelSolidFlat();
	void DrawModelSolidSmooth();
	void DrawModelSolidWithSharpness();
	void DrawModelSolidAndWireframe();
    void DrawModelHighlightOnly();
    void DrawModelTextureMapping();
    void DrawModelVertexColor();
    void DrawModelFaceColor();
	void DrawModelTransparent();
	void DrawModelFaceTexture();

    // Only render the depth of the model
    void DrawModelDepth();

    // Various auxiliary data rendering functions
    void DrawPoint();
    void DrawLine();
    void DrawPolygon();

    // Various utility rendering functions
    void DrawAxis();
    void DrawBoundingBox();
	void DrawMeshCurvature(int mode);
	void DrawVector(double scale, Coord& start, Coord& vec, Coord& normal);

	
public:
    // Constructor
    MeshModelRender();

    // Destructor
    ~MeshModelRender();

    // Initializer
    void ClearData();
    void AttachKernel(MeshModelKernel* pKernel);
    void AttachAuxData(MeshModelAuxData* pAuxData);
    
    // Get/Set functions
    int& Mode() { return m_Mode; }
    int& State() { return m_State; }
    
    Color& BbColor()  { return m_BbColor; }
    
    Color& DftVtxColor()    { return m_DftVtxColor; }
    Color& DftEdgeColor()   { return m_DftEdgeColor; }
    Color& DftFaceColor()   { return m_DftFaceColor; }
    
    Color& DftPointColor()  { return m_DftPointColor; }
    Color& DftLineColor()   { return m_DftLineColor; }
    Color& DftPolygonColor()   { return m_DftPolygonColor; }
    Color& DftElementColor()   { return m_DftElementColor; }
    
    Color& DftHltVtxColor()    { return m_DftHltVtxColor; }
    Color& DftHltEdgeColor()   { return m_DftHltEdgeColor; }
    Color& DftHltFaceColor()   { return m_DftHltFaceColor; }

    float& DftVtxSize()    { return m_DftVtxSize; }
    float& DftEdgeWidth()  { return m_DftEdgeWidth; }
    
    float& DftPointSize()  { return m_DftPointSize; }
    float& DftLineWidth()  { return m_DftLineWidth; }
    
    float& DftBdyVtxSize()     { return m_DftBdyVtxSize; }
    float& DftBdyEdgeWidth()   { return m_DftBdyEdgeWidth; }
    
    float& DftHltVtxSize()     { return m_DftHltVtxSize; }
    float& DftHltEdgeWidth()   { return m_DftHltEdgeWidth; }
    
    Material& ModelMaterial()      { return m_ModelMaterial; }
    Material& ElementMaterial()    { return m_ElementMaterial; }

    // Rendering entrance function
    void DrawModel();

	int CreateTexture(const std::string& file_name);


};