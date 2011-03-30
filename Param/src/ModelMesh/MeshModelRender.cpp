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
// MeshModelRender.cpp
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

#include "MeshModelRender.h"
#include "../Common/Utility.h"
#include <cassert>

// Constructor
MeshModelRender::MeshModelRender()
{
    kernel = NULL;
    auxdata = NULL;
}

// Destructor
MeshModelRender::~MeshModelRender()
{
}

// Initialize
void MeshModelRender::ClearData()
{
    m_Mode = RENDER_MODEL_SOLID_FLAT;

    m_State = RENDER_MODEL_AXIS | RENDER_MODEL_BOUNDING_BOX | RENDER_MODEL_LIGHTING;
    m_State |= RENDER_MODEL_POINT | RENDER_MODEL_LINE | RENDER_MODEL_POLYGON;
    
    m_BbColor = WHITE;
    
    m_DftVtxColor = RED;
    m_DftEdgeColor = GREEN;
    m_DftFaceColor = BLUE;
    
    m_DftPointColor = LIGHT_RED;
    m_DftLineColor = LIGHT_GREEN;
    m_DftPolygonColor = LIGHT_BLUE;
    m_DftElementColor = LIGHT_GREY;
    
    m_DftHltVtxColor = YELLOW;
    m_DftHltEdgeColor = CYAN;
    m_DftFaceColor = MAGENTA;
    
    m_DftVtxSize = 3.0f;
    m_DftEdgeWidth = 1.5f;
    
    m_DftPointSize = 6.0f;
    m_DftLineWidth = 3.0f;
    
    m_DftBdyVtxSize = 3.0f;
    m_DftBdyEdgeWidth = 2.0f;
    
    m_DftHltVtxSize = 6.0f;
    m_DftHltEdgeWidth = 4.0f;
        
	m_ModelMaterial.m_Ambient[0] = 0.2f;
	m_ModelMaterial.m_Ambient[1] = 0.2f;
	m_ModelMaterial.m_Ambient[2] = 0.2f;
	m_ModelMaterial.m_Ambient[3] = 1.0f;
	
	m_ModelMaterial.m_Diffuse[0] = 0.8f;
	m_ModelMaterial.m_Diffuse[1] = 0.8f;
	m_ModelMaterial.m_Diffuse[2] = 0.8f;
	m_ModelMaterial.m_Diffuse[3] = 1.0f;
	
	m_ModelMaterial.m_Specular[0] = 0.0f;
	m_ModelMaterial.m_Specular[1] = 0.0f;
	m_ModelMaterial.m_Specular[2] = 0.0f;
	m_ModelMaterial.m_Specular[3] = 1.0f;
	
	m_ModelMaterial.m_Emissive[0] = 0.0f;
	m_ModelMaterial.m_Emissive[1] = 0.0f;
	m_ModelMaterial.m_Emissive[2] = 0.0f;
	m_ModelMaterial.m_Emissive[3] = 1.0f;

	 m_ModelMaterial.m_Shininess = 0.0f;
	 m_ModelMaterial.m_Transparency = 0.0f;

    m_ElementMaterial = m_ModelMaterial;
}

// Initializer
void MeshModelRender::AttachKernel(MeshModelKernel* pKernel)
{
    assert(pKernel != NULL);
    kernel = pKernel;

    ClearData();
}

void MeshModelRender::AttachAuxData(MeshModelAuxData* pAuxData)
{
    assert(pAuxData != NULL);
    auxdata = pAuxData;
}

// Rendering entrance function
void MeshModelRender::DrawModel()
{
	// set the default material here.
	m_DefaultMaterial.SetMaterial();

    Utility util;

    // Whether show axis
    if(util.IsSetFlag(m_State, RENDER_MODEL_AXIS))
        DrawAxis();

    // Whether show bounding box
    if(util.IsSetFlag(m_State, RENDER_MODEL_BOUNDING_BOX))
        DrawBoundingBox();

    // Whether lighting
    if(util.IsSetFlag(m_State, RENDER_MODEL_LIGHTING))
        glEnable(GL_LIGHTING);
    else
        glDisable(GL_LIGHTING);

    glEnable(GL_POINT_SMOOTH);
    glEnable(GL_LINE_SMOOTH);
    glEnable(GL_POLYGON_SMOOTH);
    glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);
    glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
    glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST);
    
    // Whether per-element color
    if(util.IsSetFlag(m_State, RENDER_MODEL_COLOR))
    {
        glEnable(GL_COLOR_MATERIAL);
        glColorMaterial(GL_FRONT_AND_BACK, GL_DIFFUSE);
    }
    else
        glDisable(GL_COLOR_MATERIAL);
    
    switch(m_Mode) {
    case RENDER_MODEL_VERTICES:
        DrawModelDepth();
        DrawModelVertices();
    	break;
    case RENDER_MODEL_WIREFRAME:
        DrawModelDepth();
        DrawModelWireframe();
    	break;
    case RENDER_MODEL_SOLID_FLAT:
        DrawModelSolidFlat();
    	break;
    case RENDER_MODEL_SOLID_SMOOTH:
        DrawModelSolidSmooth();
    	break;
    case RENDER_MODEL_SOLID_WITH_SHARPNESS:
        DrawModelSolidWithSharpness();
    	break;
    case RENDER_MODEL_SOLID_AND_WIREFRAME:
        DrawModelSolidAndWireframe();
    	break;
    case RENDER_MODEL_HIGHLIGHT_ONLY:
        DrawModelHighlightOnly();
    	break;
    case RENDER_MODEL_TEXTURE_MAPPING:
        DrawModelTextureMapping();
		break;
    case RENDER_MODEL_VERTEX_COLOR:
        DrawModelVertexColor();
        break;
    case RENDER_MODEL_FACE_COLOR:
        DrawModelFaceColor();
        break;
	case RENDER_MODEL_TRANSPARENT:
		DrawModelTransparent();
		break;
	case RENDER_MODEL_PARAM_TEXTURE:
		DrawModelFaceTexture();
		break;
    }

    glDisable(GL_LIGHTING);

    if(util.IsSetFlag(m_State, RENDER_MODEL_POINT))
        DrawPoint();
    if(util.IsSetFlag(m_State, RENDER_MODEL_LINE))
        DrawLine();
    if(util.IsSetFlag(m_State, RENDER_MODEL_POLYGON))
        DrawPolygon();

	// draw kmax and kmin here.
	 if(util.IsSetFlag(m_State, RENDER_MODEL_KMAX_CURVATURE)
		 || util.IsSetFlag(m_State, RENDER_MODEL_KMIN_CURVATURE))
	 {
		 DrawMeshCurvature(m_State);
	 }

    glEnable(GL_LIGHTING);
}

// Various rendering functions
void MeshModelRender::DrawModelVertices()
{
    glDisable(GL_LIGHTING);
    glDepthFunc(GL_LEQUAL);

    Color c = RED;
    glColor3d(c[0], c[1], c[2]);
    glPointSize(m_DftVtxSize);

    CoordArray& vCoord = kernel->GetVertexInfo().GetCoord();
    FlagArray& vFlag = kernel->GetVertexInfo().GetFlag();
    size_t nVertex = vCoord.size();
    glBegin(GL_POINTS);
    for(size_t i = 0; i < nVertex; ++ i)
    {
        if(util.IsSetFlag(vFlag[i], FLAG_INVALID))
            continue;
        Coord& v = vCoord[i];
        glVertex3d(v[0], v[1], v[2]);
    }
    glEnd();
    
    glDepthFunc(GL_LESS);
    glEnable(GL_LIGHTING);
}

void MeshModelRender::DrawModelWireframe()
{
    glDisable(GL_LIGHTING);
    glDepthFunc(GL_LEQUAL);
    
    glEnable(GL_LINE_SMOOTH);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
    glLineWidth(m_DftEdgeWidth);
    
    CoordArray& vCoord = kernel->GetVertexInfo().GetCoord();
    PolyIndexArray& fIndex = kernel->GetFaceInfo().GetIndex();
    NormalArray& fNormal = kernel->GetFaceInfo().GetNormal();
    
    size_t nFace = fIndex.size();
    size_t i, j, n;
    
    Color c = BLUE*0.70;
    glColor3d(c[0], c[1], c[2]);
    
    glBegin(GL_LINES);
    for(i = 0; i < nFace; ++ i)
    {
        IndexArray& f = fIndex[i];
        n = f.size();
        for(j = 0; j < n; ++ j)
        {
            Coord& v1 = vCoord[f[j]];
            Coord& v2 = vCoord[f[(j+1)%n]];
            glVertex3d(v1[0], v1[1], v1[2]);
            glVertex3d(v2[0], v2[1], v2[2]);
        }
    }
    glEnd();
    
    glHint(GL_LINE_SMOOTH_HINT, GL_DONT_CARE);
    glDisable(GL_BLEND);
    glDisable(GL_LINE_SMOOTH);
    glEnable(GL_LIGHTING);
    glDepthFunc(GL_LESS);
    glEnable(GL_LIGHTING);
}

void MeshModelRender::DrawModelSolidFlat()
{
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glEnable(GL_POLYGON_OFFSET_FILL);
	glPolygonOffset(2.0, 2.0);
	
    CoordArray& vCoord = kernel->GetVertexInfo().GetCoord();
    PolyIndexArray& fIndex = kernel->GetFaceInfo().GetIndex();
    NormalArray& fNormal = kernel->GetFaceInfo().GetNormal();

    size_t nFace = fIndex.size();
    size_t i, j, n;

    if(!kernel->GetModelInfo().IsGeneralMesh())
    {
        if(kernel->GetModelInfo().IsTriMesh())  // Triangle mesh
            glBegin(GL_TRIANGLES);
        else if(kernel->GetModelInfo().IsQuadMesh())    // Quadangle mesh
            glBegin(GL_QUADS);
        else
            return;

        for(i = 0; i < nFace; ++ i)
        {
            IndexArray& face = fIndex[i];
            n = face.size();
            for(j = 0; j < n; ++ j)
            {
                VertexID& vID = face[j];
                Coord& v = vCoord[vID];
                Normal& n = fNormal[i];
                glNormal3d(n[0], n[1], n[2]);
                glVertex3d(v[0], v[1], v[2]);
            }
        }
        glEnd();
    }

    else    // General polygonal mesh
    {
        for(i = 0; i < nFace; ++ i)
        {
            IndexArray& face = fIndex[i];
            n = face.size();
            glBegin(GL_POLYGON);
            for(j = 0; j < n; ++ j)
            {
                VertexID vID = face[j];
                Coord& v = vCoord[vID];
                Normal& n = fNormal[i];
                glNormal3d(n[0], n[1], n[2]);
                glVertex3d(v[0], v[1], v[2]);
            }
            glEnd();
        }
    }
    glDisable(GL_POLYGON_OFFSET_FILL);
}

void MeshModelRender::DrawModelTransparent()
{
	glDisable(GL_LIGHTING);
	glColor4f(0.3f,0.3f,0.3f,0.5f);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE);  // 基于源象素alpha通道值的半透明混合函数 ( 新增 )
	glEnable(GL_BLEND);     // 打开混合
    glDisable(GL_DEPTH_TEST); // 关闭深度测试
	DrawModelSolidSmooth();
	glEnable(GL_DEPTH_TEST);
	glDisable(GL_BLEND);
	glEnable(GL_LIGHTING);
}


void MeshModelRender::DrawModelFaceTexture()
{
	glEnable(GL_POLYGON_SMOOTH);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glEnable(GL_POLYGON_OFFSET_FILL);
	glPolygonOffset(2.0, 2.0);
	glDepthFunc(GL_LEQUAL);
	//glShadeModel(GL_SMOOTH);

	CoordArray& vCoord = kernel->GetVertexInfo().GetCoord();
	TexCoordArray& vTexCoord = kernel->GetVertexInfo().GetTexCoord();
	PolyTexCoordArray& fTexCoord = kernel->GetFaceInfo().GetTexCoord();
	NormalArray& fNormal = kernel->GetFaceInfo().GetNormal();
	PolyIndexArray& fIndex = kernel->GetFaceInfo().GetIndex();
	PolyIndexArray& ftIndex = kernel->GetFaceInfo().GetTexIndex();
    
	if(fTexCoord.size() ==0) return;

	size_t nFace = fIndex.size();
	size_t i, j, n;

	// set texture env here.
	glEnable(GL_TEXTURE_2D);
	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
	//glBindTexture(GL_TEXTURE_1D, texName);

	if(!kernel->GetModelInfo().IsGeneralMesh())
	{
		if(kernel->GetModelInfo().IsTriMesh())  // Triangle mesh
			glBegin(GL_TRIANGLES);
		else if(kernel->GetModelInfo().IsQuadMesh())    // Quadangle mesh
			glBegin(GL_QUADS);
		else
			return;

		for(i = 0; i < nFace; ++ i)
		{
			IndexArray& face = fIndex[i];
            //			IndexArray& tex_index = ftIndex[i];
			TexCoordArray& f_tex = fTexCoord[i];
			Normal& fn = fNormal[i];

			n = face.size();
			for(j = 0; j < n; ++ j)
			{
				VertexID& vID = face[j];
				Coord& v = vCoord[vID];
                //	const TexCoord& tex = vTexCoord[tex_index[j]];
				glNormal3d(fn[0], fn[1], fn[2]);
				//glTexCoord2d(tex[0], tex[1]);
				glTexCoord2d(f_tex[j][0], f_tex[j][1]);
				glVertex3d(v[0], v[1], v[2]);
			}
		}
		glEnd();
	}

	glDisable(GL_TEXTURE_2D);
	glDisable(GL_POLYGON_OFFSET_FILL);
	glDisable(GL_POLYGON_SMOOTH);
}

void MeshModelRender::DrawModelSolidSmooth()
{
    glEnable(GL_POLYGON_SMOOTH);
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glEnable(GL_POLYGON_OFFSET_FILL);
	glPolygonOffset(2.0, 2.0);
	
    CoordArray& vCoord = kernel->GetVertexInfo().GetCoord();
    NormalArray& vNormal = kernel->GetVertexInfo().GetNormal();
    PolyIndexArray& fIndex = kernel->GetFaceInfo().GetIndex();

    size_t nFace = fIndex.size();
    size_t i, j, n;

    if(!kernel->GetModelInfo().IsGeneralMesh())
    {
        if(kernel->GetModelInfo().IsTriMesh())  // Triangle mesh
            glBegin(GL_TRIANGLES);
        else if(kernel->GetModelInfo().IsQuadMesh())    // Quadangle mesh
            glBegin(GL_QUADS);
        else
            return;

        for(i = 0; i < nFace; ++ i)
        {
            IndexArray& face = fIndex[i];
            n = face.size();
            for(j = 0; j < n; ++ j)
            {
                VertexID& vID = face[j];
                Coord& v = vCoord[vID];
                Normal& n = vNormal[vID];
                glNormal3d(n[0], n[1], n[2]);
                glVertex3d(v[0], v[1], v[2]);
            }
        }
        glEnd();
    }

    else    // General polygonal mesh
    {
        for(i = 0; i < nFace; ++ i)
        {
            IndexArray& face = fIndex[i];
            n = face.size();
            glBegin(GL_POLYGON);
            for(j = 0; j < n; ++ j)
            {
                VertexID vID = face[j];
                Coord& v = vCoord[vID];
                Normal& n = vNormal[vID];
                glNormal3d(n[0], n[1], n[2]);
                glVertex3d(v[0], v[1], v[2]);
            }
            glEnd();
        }
    }
    glDisable(GL_POLYGON_OFFSET_FILL);
    glDisable(GL_POLYGON_SMOOTH);
}

void MeshModelRender::DrawModelSolidWithSharpness()
{
}

void MeshModelRender::DrawModelSolidAndWireframe()
{
    // Draw solid smooth model
    glEnable(GL_POLYGON_SMOOTH);
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glEnable(GL_POLYGON_OFFSET_FILL);
	glPolygonOffset(2.0, 2.0);
	
    CoordArray& vCoord = kernel->GetVertexInfo().GetCoord();
    NormalArray& vNormal = kernel->GetVertexInfo().GetNormal();
    PolyIndexArray& fIndex = kernel->GetFaceInfo().GetIndex();

    size_t nFace = fIndex.size();
    size_t i, j, n;

    if(!kernel->GetModelInfo().IsGeneralMesh())
    {
        if(kernel->GetModelInfo().IsTriMesh())  // Triangle mesh
            glBegin(GL_TRIANGLES);
        else if(kernel->GetModelInfo().IsQuadMesh())    // Quadangle mesh
            glBegin(GL_QUADS);
        else
            return;

        for(i = 0; i < nFace; ++ i)
        {
            IndexArray& face = fIndex[i];
            n = face.size();
            for(j = 0; j < n; ++ j)
            {
                VertexID& vID = face[j];
                Coord& v = vCoord[vID];
                Normal& n = vNormal[vID];
                glNormal3d(n[0], n[1], n[2]);
                glVertex3d(v[0], v[1], v[2]);
            }
        }
        glEnd();
    }

    else    // General polygonal mesh
    {
        for(i = 0; i < nFace; ++ i)
        {
            IndexArray& face = fIndex[i];
            n = face.size();
            glBegin(GL_POLYGON);
            for(j = 0; j < n; ++ j)
            {
                VertexID vID = face[j];
                Coord& v = vCoord[vID];
                Normal& n = vNormal[vID];
                glNormal3d(n[0], n[1], n[2]);
                glVertex3d(v[0], v[1], v[2]);
            }
            glEnd();
        }
    }
    glDisable(GL_POLYGON_OFFSET_FILL);
    glDisable(GL_POLYGON_SMOOTH);

    // Draw wireframe
    glDisable(GL_LIGHTING);
    glDepthFunc(GL_LEQUAL);

    glEnable(GL_LINE_SMOOTH);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
    glLineWidth(m_DftEdgeWidth);
    
    Color c = BLUE*0.70;
    glColor3d(c[0], c[1], c[2]);

    glBegin(GL_LINES);
    for(i = 0; i < nFace; ++ i)
    {
        IndexArray& f = fIndex[i];
        n = f.size();
        for(j = 0; j < n; ++ j)
        {
            Coord& v1 = vCoord[f[j]];
            Coord& v2 = vCoord[f[(j+1)%n]];
            glVertex3d(v1[0], v1[1], v1[2]);
            glVertex3d(v2[0], v2[1], v2[2]);
        }
    }
    glEnd();

    glHint(GL_LINE_SMOOTH_HINT, GL_DONT_CARE);
    glDisable(GL_BLEND);
    glDisable(GL_LINE_SMOOTH);
    glEnable(GL_LIGHTING);
    glDepthFunc(GL_LESS);
    glEnable(GL_LIGHTING);
}

void MeshModelRender::DrawModelHighlightOnly()
{
}

void MeshModelRender::DrawModelTextureMapping()
{
	CoordArray& vCoord = kernel->GetVertexInfo().GetCoord();
	TexCoordArray& tCoord = kernel->GetVertexInfo().GetTexCoord();
	NormalArray& vNormal = kernel->GetVertexInfo().GetNormal();
	PolyIndexArray& fIndex = kernel->GetFaceInfo().GetIndex();

	if (tCoord.empty())
	{
		return;
	}

	glEnable(GL_POLYGON_SMOOTH);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glEnable(GL_POLYGON_OFFSET_FILL);
	glPolygonOffset(2.0, 2.0);
	glDepthFunc(GL_LEQUAL);
	//glShadeModel(GL_SMOOTH);
	size_t nFace = fIndex.size();
	size_t i, j, n;

	// set texture env here.
	glEnable(GL_TEXTURE_2D);
	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
	//glBindTexture(GL_TEXTURE_1D, texName);

	if(!kernel->GetModelInfo().IsGeneralMesh())
	{
		if(kernel->GetModelInfo().IsTriMesh())  // Triangle mesh
			glBegin(GL_TRIANGLES);
		else if(kernel->GetModelInfo().IsQuadMesh())    // Quadangle mesh
			glBegin(GL_QUADS);
		else
			return;

		for(i = 0; i < nFace; ++ i)
		{
			IndexArray& face = fIndex[i];
			n = face.size();
			for(j = 0; j < n; ++ j)
			{
				VertexID& vID = face[j];
				Coord& v = vCoord[vID];
				Coord2D& t = tCoord[vID];
				Normal& n = vNormal[vID];
				glNormal3d(n[0], n[1], n[2]);
				glTexCoord2d(t[0], t[1]);
				glVertex3d(v[0], v[1], v[2]);
			}
		}
		glEnd();
	}
	else    // General polygonal mesh
	{
		for(i = 0; i < nFace; ++ i)
		{
			IndexArray& face = fIndex[i];
			n = face.size();
			glBegin(GL_POLYGON);
			for(j = 0; j < n; ++ j)
			{
				VertexID vID = face[j];
				Coord& v = vCoord[vID];
				Coord2D& t = tCoord[vID];
				Normal& n = vNormal[vID];
				glNormal3d(n[0], n[1], n[2]);
				glTexCoord2d(t[0], t[1]);
				glVertex3d(v[0], v[1], v[2]);
			}
			glEnd();
		}
	}

	glDisable(GL_TEXTURE_2D);
	glDisable(GL_POLYGON_OFFSET_FILL);
	glDisable(GL_POLYGON_SMOOTH);
}

void MeshModelRender::DrawModelVertexColor()
{
    CoordArray& vCoord = kernel->GetVertexInfo().GetCoord();
    NormalArray& vNormal = kernel->GetVertexInfo().GetNormal();
    ColorArray& vColor = kernel->GetVertexInfo().GetColor();
    PolyIndexArray& fIndex = kernel->GetFaceInfo().GetIndex();

    if(vColor.size() != vCoord.size())
        return;

    // Get previous material diffuse components
    float ambient[4], diffuse[4], specular[4], emission[4], shininess;
    glGetMaterialfv(GL_FRONT, GL_AMBIENT, ambient);
    glGetMaterialfv(GL_FRONT, GL_DIFFUSE, diffuse);
    glGetMaterialfv(GL_FRONT, GL_SPECULAR, specular);
    glGetMaterialfv(GL_FRONT, GL_EMISSION, emission);
    glGetMaterialfv(GL_FRONT, GL_SHININESS, &shininess);

	glDisable(GL_LIGHTING);
	glShadeModel(GL_SMOOTH);

    //glEnable(GL_POLYGON_SMOOTH);
    //glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	//glEnable(GL_POLYGON_OFFSET_FILL);
	//glPolygonOffset(2.0, 2.0);
    glEnable(GL_COLOR_MATERIAL);
    glColorMaterial(GL_FRONT_AND_BACK, GL_DIFFUSE);

    size_t nFace = fIndex.size();
    size_t i, j, n;

    if(!kernel->GetModelInfo().IsGeneralMesh())
    {
        if(kernel->GetModelInfo().IsTriMesh())  // Triangle mesh
            glBegin(GL_TRIANGLES);
        else if(kernel->GetModelInfo().IsQuadMesh())    // Quadangle mesh
            glBegin(GL_QUADS);
        else
            return;

        for(i = 0; i < nFace; ++ i)
        {
            IndexArray& face = fIndex[i];
            n = face.size();
            for(j = 0; j < n; ++ j)
            {
                VertexID& vID = face[j];
                Coord& v = vCoord[vID];
                Normal& n = vNormal[vID];
                Color& c = vColor[vID];
                glColor3d (c[0], c[1], c[2]);
                glNormal3d(n[0], n[1], n[2]);
                glVertex3d(v[0], v[1], v[2]);
            }
        }
        glEnd();
    }

    else    // General polygonal mesh
    {
        for(i = 0; i < nFace; ++ i)
        {
            IndexArray& face = fIndex[i];
            n = face.size();
            glBegin(GL_POLYGON);
            for(j = 0; j < n; ++ j)
            {
                VertexID vID = face[j];
                Coord& v = vCoord[vID];
                Normal& n = vNormal[vID];
                Color& c = vColor[vID];
                glColor3d (c[0], c[1], c[2]);
                glNormal3d(n[0], n[1], n[2]);
                glVertex3d(v[0], v[1], v[2]);
            }
            glEnd();
        }
    }

    glDisable(GL_COLOR_MATERIAL);
    glDisable(GL_POLYGON_OFFSET_FILL);
    glDisable(GL_POLYGON_SMOOTH);
	glEnable(GL_LIGHTING);

    // Reset previous material diffuse components
    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, ambient);
    glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, diffuse);
    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, specular);
    glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, emission);
    glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, shininess);
}

void MeshModelRender::DrawModelFaceColor()
{
    CoordArray& vCoord = kernel->GetVertexInfo().GetCoord();
    NormalArray& vNormal = kernel->GetVertexInfo().GetNormal();
    ColorArray& fColor = kernel->GetFaceInfo().GetColor();
    PolyIndexArray& fIndex = kernel->GetFaceInfo().GetIndex();

    if(fColor.size() != fIndex.size())
        return;

    // Get previous material diffuse components
    float ambient[4], diffuse[4], specular[4], emission[4], shininess;
    glGetMaterialfv(GL_FRONT, GL_AMBIENT, ambient);
    glGetMaterialfv(GL_FRONT, GL_DIFFUSE, diffuse);
    glGetMaterialfv(GL_FRONT, GL_SPECULAR, specular);
    glGetMaterialfv(GL_FRONT, GL_EMISSION, emission);
    glGetMaterialfv(GL_FRONT, GL_SHININESS, &shininess);


	glDisable(GL_LIGHTING);
	glShadeModel(GL_SMOOTH);

//  glEnable(GL_POLYGON_SMOOTH);
//  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
// 	glEnable(GL_POLYGON_OFFSET_FILL);
// 	glPolygonOffset(2.0, 2.0);
    glEnable(GL_COLOR_MATERIAL);
    glColorMaterial(GL_FRONT, GL_DIFFUSE);

    size_t nFace = fIndex.size();
    size_t i, j, n;

    if(!kernel->GetModelInfo().IsGeneralMesh())
    {
        if(kernel->GetModelInfo().IsTriMesh())  // Triangle mesh
            glBegin(GL_TRIANGLES);
        else if(kernel->GetModelInfo().IsQuadMesh())    // Quadangle mesh
            glBegin(GL_QUADS);
        else
            return;

        for(i = 0; i < nFace; ++ i)
        {			
			IndexArray& face = fIndex[i];
            n = face.size();
            Color& c = fColor[i];
            glColor3d(c[0], c[1], c[2]);			
            for(j = 0; j < n; ++ j)
            {
                VertexID& vID = face[j];
                Coord& v = vCoord[vID];
                Normal& n = vNormal[vID];
                glNormal3d(n[0], n[1], n[2]);
                glVertex3d(v[0], v[1], v[2]);
            }
        }
        glEnd();
    }

    else    // General polygonal mesh
    {
        for(i = 0; i < nFace; ++ i)
        {
            IndexArray& face = fIndex[i];
            n = face.size();
            Color& c = fColor[i];
            glColor3d(c[0], c[1], c[2]);
            glBegin(GL_POLYGON);
            for(j = 0; j < n; ++ j)
            {
                VertexID vID = face[j];
                Coord& v = vCoord[vID];
                Normal& n = vNormal[vID];
                glNormal3d(n[0], n[1], n[2]);
                glVertex3d(v[0], v[1], v[2]);
            }
            glEnd();
        }
    }
    glDisable(GL_COLOR_MATERIAL);
    glDisable(GL_POLYGON_OFFSET_FILL);
    glDisable(GL_POLYGON_SMOOTH);
	glEnable(GL_LIGHTING);

    // Reset previous material diffuse components
    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, ambient);
    glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, diffuse);
    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, specular);
    glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, emission);
    glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, shininess);
}

// Various utility rendering functions
void MeshModelRender::DrawAxis()
{
}

void MeshModelRender::DrawBoundingBox()
{
	Coord BoxMin, BoxMax, BoxDim, Center;
    
    kernel->GetModelInfo().GetBoundingBox(BoxMin, BoxMax, BoxDim);
	
	// Draw bouding box
    glDisable(GL_LIGHTING);
	
    glEnable(GL_LINE_SMOOTH);
    glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
    glLineWidth(m_DftEdgeWidth);
    
    Color c = GREEN*0.8;
    glColor3d(c[0], c[1], c[2]);
	
    glBegin(GL_LINES);
    
	// 
    glVertex3d(BoxMin[0], BoxMin[1], BoxMax[2]);
	glVertex3d(BoxMax[0], BoxMin[1], BoxMax[2]);
	
	glVertex3d(BoxMax[0], BoxMin[1], BoxMax[2]);
	glVertex3d(BoxMax[0], BoxMin[1], BoxMin[2]);
	
	glVertex3d(BoxMax[0], BoxMin[1], BoxMin[2]);
	glVertex3d(BoxMin[0], BoxMin[1], BoxMin[2]);
	
	glVertex3d(BoxMin[0], BoxMin[1], BoxMin[2]);
	glVertex3d(BoxMin[0], BoxMin[1], BoxMax[2]);
	
	//
	glVertex3d(BoxMin[0], BoxMax[1], BoxMax[2]);
	glVertex3d(BoxMax[0], BoxMax[1], BoxMax[2]);
	
	glVertex3d(BoxMax[0], BoxMax[1], BoxMax[2]);
	glVertex3d(BoxMax[0], BoxMax[1], BoxMin[2]);
	
	glVertex3d(BoxMax[0], BoxMax[1], BoxMin[2]);
	glVertex3d(BoxMin[0], BoxMax[1], BoxMin[2]);
	
	glVertex3d(BoxMin[0], BoxMax[1], BoxMin[2]);
	glVertex3d(BoxMin[0], BoxMax[1], BoxMax[2]);
	
	//
	glVertex3d(BoxMin[0], BoxMin[1], BoxMax[2]);
	glVertex3d(BoxMin[0], BoxMax[1], BoxMax[2]);
	
	glVertex3d(BoxMax[0], BoxMin[1], BoxMax[2]);
	glVertex3d(BoxMax[0], BoxMax[1], BoxMax[2]);
	
    glVertex3d(BoxMax[0], BoxMin[1], BoxMin[2]);
	glVertex3d(BoxMax[0], BoxMax[1], BoxMin[2]);
	
	glVertex3d(BoxMin[0], BoxMin[1], BoxMin[2]);
	glVertex3d(BoxMin[0], BoxMax[1], BoxMin[2]);
	
    glEnd();
	
    glHint(GL_LINE_SMOOTH_HINT, GL_DONT_CARE);
    glDisable(GL_LINE_SMOOTH);
    glEnable(GL_LIGHTING);
}

// Various auxiliary data rendering functions
void MeshModelRender::DrawPoint()
{
    glEnable(GL_POINT_SMOOTH);
    glDepthFunc(GL_LEQUAL);
    
    glPointSize(m_DftPointSize);

    PointInfo& pInfo = auxdata->GetPointInfo();
    CoordArray& pCoord = pInfo.GetCoord();
    ColorArray& pColor = pInfo.GetColor();

    size_t n = pCoord.size();
    glBegin(GL_POINTS);
    for(size_t i = 0; i < n; ++ i)
    {
        Coord& p = pCoord[i];
        Color& c = pColor[i];
        glColor3d(c[0], c[1], c[2]);
        glVertex3d(p[0], p[1], p[2]);
    }
    glEnd();

    glDepthFunc(GL_LESS);
    glDisable(GL_POINT_SMOOTH);
}

void MeshModelRender::DrawLine()
{
    glEnable(GL_LINE_SMOOTH);
    glDepthFunc(GL_LEQUAL);

    glLineWidth(m_DftLineWidth);

    LineInfo& lInfo = auxdata->GetLineInfo();
    CoordArray& lfCoord = lInfo.GetfCoord();
    CoordArray& lsCoord = lInfo.GetsCoord();
    ColorArray& lColor = lInfo.GetColor();

    size_t n = lfCoord.size();
    glBegin(GL_LINES);
    for(size_t i = 0; i < n ; ++ i)
    {
        Coord& fp = lfCoord[i];
        Coord& sp = lsCoord[i];
        Color& c = lColor[i];
        glColor3d(c[0], c[1], c[2]);
        glVertex3d(fp[0], fp[1], fp[2]);
        glVertex3d(sp[0], sp[1], sp[2]);
    }
    glEnd();

    glDepthFunc(GL_LESS);
    glDisable(GL_LINE_SMOOTH);
}

void MeshModelRender::DrawPolygon()
{
    glEnable(GL_POLYGON_SMOOTH);
    glEnable(GL_POLYGON_OFFSET_FILL);

    glPolygonOffset(1.0, 1.0);

    PolygonInfo& pInfo = auxdata->GetPolygonInfo();
    PolyCoordArray& pCoords = pInfo.GetCoords();
    ColorArray& pColor = pInfo.GetColor();

    size_t n = pCoords.size();
    for(size_t i = 0; i < n; ++ i)
    {
        Color& c = pColor[i];
        glColor3d(c[0], c[1], c[2]);

        CoordArray& pv = pCoords[i];
        size_t m = pv.size();
        glBegin(GL_POLYGON);
        for(size_t j = 0; j < m; ++ j)
        {
            Coord& v = pv[j];
            glVertex3d(v[0], v[1], v[2]);
        }
        glEnd();
    }

    glDisable(GL_POLYGON_OFFSET_FILL);
    glDisable(GL_POLYGON_SMOOTH);
}

// Rendering object depth for clipping
void MeshModelRender::DrawModelDepth()
{
    glColorMask(GL_FALSE, GL_FALSE, GL_FALSE, GL_FALSE);
    
    DrawModelSolidFlat();

    glColorMask(GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);
}
void MeshModelRender::DrawMeshCurvature(int mode)
{
	double scale = 0.75 * kernel->GetModelInfo().GetAvgEdgeLength();

	CoordArray& vCoord = kernel->GetVertexInfo().GetCoord();
	CoordArray& vNormal = kernel->GetVertexInfo().GetNormal();
	CurvatureArray& vCurvature = kernel->GetVertexInfo().GetCurvatures();
	
	glDisable(GL_LIGHTING);
	glEnable(GL_POLYGON_OFFSET_LINE);
	glPolygonOffset(3.0f, -1.0f);
	glEnable(GL_LINE_SMOOTH);
	
	size_t nVertex = vCoord.size();
    size_t j;
	
	for (j = 0; j < nVertex; j++)
	{
		Coord& start = vCoord[j];
		Coord& normal = vNormal[j];
		CCurvature& curv = vCurvature[j];
		
		// draw the kmax curvature here
		if (util.IsSetFlag(mode, RENDER_MODEL_KMAX_CURVATURE))
		{
			glColor3d(1.0f,0.0f,1.0f);
			
			DrawVector(scale, start, curv.m_direction_kmax, normal);
		}
		
		// draw the kmin curvature here.
		if (util.IsSetFlag(mode, RENDER_MODEL_KMIN_CURVATURE))
		{
			glColor3d(0.0f,0.0f,1.0f);
			
			DrawVector(scale, start, curv.m_direction_kmin, normal);
		}
	}
	
	glDisable(GL_LINE_SMOOTH);
	glDisable(GL_POLYGON_OFFSET_LINE);
	glEnable(GL_LIGHTING);
}
void MeshModelRender::DrawVector(double scale, Coord& start, Coord& vec, Coord& normal)
{
	double x1 = (double) ( scale * vec[0]);
	double y1 = (double) ( scale * vec[1]);
	double z1 = (double) ( scale * vec[2]);
	
	double x2 = x1 * 0.8;
	double y2 = y1 * 0.8;
	double z2 = z1 * 0.8;
	
	Coord isoV = cross(normal, vec).unit();
    Coord p0 = Coord(x2, y2, z2) + isoV * scale * 0.15;
	Coord p1 = Coord(x2, y2, z2) - isoV * scale * 0.15;
	
	glPushMatrix();
	glTranslated(start[0], start[1], start[2]);
	
	glLineWidth(1.0f);
	glBegin(GL_LINES);
	
	glVertex3d(0,0,0);
	glVertex3d(x1,y1,z1);
	
	glEnd();
	
	glBegin(GL_LINES);
	
	
	glVertex3d(p0[0], p0[1], p0[2]);
	glVertex3d(x1,y1,z1);
	
	glEnd();
	
	glBegin(GL_LINES);
	
	
	glVertex3d(p1[0], p1[1], p1[2]);
	glVertex3d(x1,y1,z1);
	
	glEnd();
	
	glPopMatrix();                              
}


int MeshModelRender::CreateTexture(const std::string& texture_file_name)
{

#ifdef WIN32
	GLuint texName;

	AUX_RGBImageRec* texture_image = auxDIBImageLoadA((LPCSTR)texture_file_name.c_str());

	glGenTextures(1, &texName);					

	glBindTexture(GL_TEXTURE_2D, texName);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_NEAREST);
	glTexImage2D(GL_TEXTURE_2D, 0, 3, texture_image->sizeX, 
	texture_image->sizeY, 0, GL_RGB, GL_UNSIGNED_BYTE, 
	texture_image->data);

	if (texture_image)
	{
	 	if (texture_image->data)	free(texture_image->data);			
	 	free(texture_image);					
	}
#endif
	return 0;
}
