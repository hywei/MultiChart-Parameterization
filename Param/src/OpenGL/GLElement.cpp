
#include "GLElement.h"
#include <math.h>
#include <assert.h>
#include <iostream>
#include <QtOpenGL/QGLWidget>
#include <cassert>
#include <fstream>


/* ================== Material Information ================== */

// Constructor
GLMaterial::GLMaterial()
{
	ClearData();
}

GLMaterial::GLMaterial(const GLMaterial& material)
{
	int i;
	
	for(i = 0; i < 4; ++ i)
	{
        this->ambient[i] = material.ambient[i];
	}
	
	for(i = 0; i < 4; ++ i)
	{
        this->diffuse[i] = material.diffuse[i];
	}
	
	for(i = 0; i < 4; ++ i)
	{
        this->specular[i] = material.specular[i];
	}
	
	for(i = 0; i < 4; ++ i)
	{
        this->emission[i] = material.emission[i];
	}
	
	this->shininess = material.shininess;
	this->face = material.face;
}

// Destructor
GLMaterial::~GLMaterial()
{
}


// Initializer
void GLMaterial::ClearData()
{
    ambient[0] = 0.2f;
    ambient[1] = 0.2f;
    ambient[2] = 0.2f;
    ambient[3] = 1.0f;

    diffuse[0] = 0.8f;
    diffuse[1] = 0.738f;
    diffuse[2] = 0.738f;
    diffuse[3] = 1.0f;

    specular[0] = 0.05f;
    specular[1] = 0.0f;
    specular[2] = 0.0f;
    specular[3] = 1.0f;

    emission[0] = 0.0f;
    emission[1] = 0.0f;
    emission[2] = 0.0f;
    emission[3] = 1.0f;

    shininess = 10.0f;

    face = GL_FRONT_AND_BACK;
}

// Set material value
void GLMaterial::SetValue(int pname, float* param)
{
    int i;
    switch(pname) {
    case MATERIAL_AMBIENT:
        for(i = 0; i < 4; ++ i)
        ambient[i] = param[i];
    	break;
    case MATERIAL_DIFFUSE:
        for(i = 0; i < 4; ++ i)
        diffuse[i] = param[i];
    	break;
    case MATERIAL_SPECULAR:
        for(i = 0; i < 4; ++ i)
        specular[i] = param[i];
    	break;
    case MATERIAL_EMISSION:
        for(i = 0; i < 4; ++ i)
        emission[i] = param[i];
    	break;
    case MATERIAL_SHININESS:
        shininess = *param;
    	break;
    }
}

// Set material
void GLMaterial::SetMaterial()
{
    glMaterialfv(face, GL_AMBIENT, ambient);
    glMaterialfv(face, GL_DIFFUSE, diffuse);
    glMaterialfv(face, GL_SPECULAR, specular);
    glMaterialfv(face, GL_EMISSION, emission);
    glMaterialf( face, GL_SHININESS, shininess);
}



/* ================== Light Information ================== */

// Constructor
GLLight::GLLight()
{
    ClearData();
}


// Destructor
GLLight::~GLLight()
{
}


// Initializer
void GLLight::ClearData()
{
    type = LIGHT_TYPE_DIRECTIONAL;
    

    ambient[0] = 0.05f;
    ambient[1] = 0.0f;
    ambient[2] = 0.0f;
    ambient[3] = 1.0f;

	
    diffuse[0] = 1.0f;
    diffuse[1] = 1.0f;
    diffuse[2] = 1.0f;
    diffuse[3] = 1.0f;

    specular[0] = 1.0f;
    specular[1] = 1.0f;
    specular[2] = 1.0f;
    specular[3] = 1.0f;

    position[0] = 0.0f;
    position[1] = 0.0f;
    position[2] = -1.0f;
    position[3] = 0.0f;

    direction[0] = 0.0f;
    direction[1] = 0.0f;
    direction[2] = 1.0f;

    exponent = 0.0f;

    cutoff = 180.0f;
    
    constant_attenuation = 1.0f;
    linear_attenuation = 0.0f;
    quadratic_attenuation = 0.0f;

    lightID = GL_LIGHT0;
	
	glEnable(lightID);
}

// Set material value
void GLLight::SetValue(int pname, float* param)
{
    int i;
    switch(pname) {
    case LIGHT_AMBIENT:
        for(i = 0; i < 4; ++ i)
        ambient[i] = param[i];
    	break;
    case LIGHT_DIFFUSE:
        for(i = 0; i < 4; ++ i)
        diffuse[i] = param[i];
    	break;
    case LIGHT_SPECULAR:
        for(i = 0; i < 4; ++ i)
        specular[i] = param[i];
    	break;
    case LIGHT_POSITION:
        for(i = 0; i < 4; ++ i)
        position[i] = param[i];
    	break;
    case LIGHT_SPOT_DIRECTION:
        for(i = 0; i < 3; ++ i)
        direction[i] = param[i];
    	break;
    case LIGHT_SPOT_EXPONENT:
        exponent = *param;
    	break;
    case LIGHT_SPOT_CUTOFF:
        cutoff = *param;
    	break;
    case LIGHT_CONSTANT_ATTENUATION:
        constant_attenuation = *param;
    	break;
    case LIGHT_LINEAR_ATTENUATION:        
        linear_attenuation = *param;
    	break;
    case LIGHT_QUADRATIC_ATTENUATION:
        quadratic_attenuation = *param;
    	break;
    }
}

// Set material
void GLLight::SetLight()
{
    glLightfv(lightID, GL_AMBIENT, ambient);
    glLightfv(lightID, GL_DIFFUSE, diffuse);
    glLightfv(lightID, GL_SPECULAR, specular);
    glLightfv(lightID, GL_POSITION, position);

    if(position[3] == 0.0f)  // Directional light
    {
        type = LIGHT_TYPE_DIRECTIONAL;

        glLightf(lightID, GL_CONSTANT_ATTENUATION, 1.0f);
        glLightf(lightID, GL_LINEAR_ATTENUATION, 0.0f);
        glLightf(lightID, GL_QUADRATIC_ATTENUATION, 0.0f);
    }
    else 
    {
        if(cutoff == 180.0f)
        {
            type = LIGHT_TYPE_POINT;
        }
        else
        {
            type = LIGHT_TYPE_SPOT;
            glLightfv(lightID, GL_SPOT_DIRECTION, direction);
            glLightf(lightID, GL_SPOT_EXPONENT, exponent);
            glLightf(lightID, GL_SPOT_CUTOFF, cutoff);
        }
        glLightf(lightID, GL_CONSTANT_ATTENUATION, constant_attenuation);
        glLightf(lightID, GL_LINEAR_ATTENUATION, linear_attenuation);
        glLightf(lightID, GL_QUADRATIC_ATTENUATION, quadratic_attenuation);
    }

    glEnable(lightID);  // Enable lightID
}
                 
              
              
/* ================== Projection Information ================== */

// Constructor
GLProjection::GLProjection()
{
    ClearData();
}

// Destructor
GLProjection::~GLProjection()
{
}

// Initializer
void GLProjection::ClearData()
{
    vv_left = vv_right = vv_top = vv_bottom = 0.0;
    fov = 40.0f;
    aspect = 1.0f;

    vv_near = 0.5;
    vv_far = 500.0;

    type = PROJECTION_TYPE_ORTHOGONAL;
//	type = PROJECTION_TYPE_PERSPECTIVE;

    vp_left = vp_right = vp_top = vp_bottom = 0;

    object_center[0] = 0.0;
    object_center[1] = 0.0;
    object_center[2] = 0.0;
    object_radius = 1.0;
    
    object_ratio = 1.3;
}

// Set object information
void GLProjection::SetObjectInfo(double cx, double cy, double cz, double r)
{
    object_center[0] = cx;
    object_center[1] = cy;
    object_center[2] = cz;
    object_radius    = r;

	vv_far = 4*r;
}

// Set projection
void GLProjection::SetProjection(int width, int height)
{
	if(height == 0) height = width+1;

    aspect = (double)width/(double)height;

    glMatrixMode( GL_PROJECTION );
	glLoadIdentity();
    if(type == PROJECTION_TYPE_ORTHOGONAL)
    {
        vv_top = object_radius*object_ratio;
        vv_bottom = -vv_top;
        vv_right = aspect*vv_top;
        vv_left = -vv_right;
        glOrtho(vv_left, vv_right, vv_bottom, vv_top, -vv_far, vv_far);
    }
    else
    {
        vv_top = vv_near*tan(fov/2.0);
        vv_bottom = -vv_top;
        vv_right = vv_top*aspect;
        vv_left = -vv_right;
        glFrustum(vv_left, vv_right, vv_bottom, vv_top, vv_near, vv_far);
    }
    
    double z_dist = (object_radius*object_ratio)/tan(fov/2.0);
  //  gluLookAt(object_center[0], object_center[1], object_center[2]+z_dist, object_center[0], object_center[1], object_center[2], 0, 1, 0);

    glViewport(0, 0, width, height);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
}

void GLProjection::SetProjection(int left, int right, int top, int bottom)
{
	std::cout<<left<<' '<<right<<' '<<top<<' '<<bottom<<std::endl;
	std::cout<<"Warning : This function does nothing!\n";
}

void GLProjection::Zoom(double scale)
{
    if(type == PROJECTION_TYPE_ORTHOGONAL)
        object_ratio *= scale;
    else
        fov *= scale;
}

/* ================== OpenGL Element Information ================== */

// Constructor
COpenGL::COpenGL()
{
    ClearData();
}

// Destructor
COpenGL::~COpenGL()
{
}

// Initializer
void COpenGL::ClearData()
{

	m_BkColor[0] = 0.0;
	m_BkColor[1] = 0.0;
	m_BkColor[2] = 0.0;
	m_BkColor[3] = 1.0;

	
    m_BkColor[0] = 0.00;
    m_BkColor[1] = 0.00;
    m_BkColor[2] = 0.5; //0.25
    m_BkColor[3] = 1.0;
	

    m_Light.ClearData();
	m_Light.SetLight();

    m_Projection.ClearData();
}

// Debug
bool COpenGL::DetectOpenGLError()
{
	GLenum errCode;
	const GLubyte *errString;

	if((errCode = glGetError()) != GL_NO_ERROR)
	{
		errString = gluErrorString(errCode);
		printf("OpenGL Error: %s\n",errString);
        assert(0);
	}
	
    return (errCode == GL_NO_ERROR);
}

void COpenGL::OnCreate(QGLWidget *glWidget)
{
	m_GLContext = glWidget;
}

void COpenGL::SetObjectInfo(double cx, double cy, double cz, double r)
{
	m_GLContext->makeCurrent();
	m_Projection.SetObjectInfo(cx, cy, cz, r);
	m_GLContext->doneCurrent();
}

void COpenGL::OnResize(int width, int height)
{
    if ((width <= 0) || (height <= 0))
        return;
	
	m_GLContext->makeCurrent();
	DetectOpenGLError();

    m_Projection.SetProjection(width, height);

    DetectOpenGLError();
    
	m_GLContext->doneCurrent();
}

void COpenGL::OnReset(int width, int height)
{
    if ((width <= 0) || (height <= 0))
        return;

	m_GLContext->makeCurrent();

	DetectOpenGLError();

    m_Projection.ClearData();
    m_Projection.SetProjection(width, height);

    DetectOpenGLError();
    
	m_GLContext->doneCurrent();
}

void COpenGL::OnInit()
{
	m_GLContext->makeCurrent();

    // Set background color
	glClearColor((float)m_BkColor[0], (float)m_BkColor[1], (float)m_BkColor[2], (float)m_BkColor[3]);
	
    DetectOpenGLError();

	// Set lights
//     glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, 1);
//     m_Light.SetLight();
// 	   glEnable(GL_LIGHTING);

//    DetectOpenGLError();

    // Set material
//    m_Material.SetMaterial();

//    DetectOpenGLError();

    // Enable depth test
	glEnable(GL_DEPTH_TEST);

	m_GLContext->doneCurrent();

	m_ModelViewMatrix.MakeUnitMatrix(4);
}

void COpenGL::OnBeginPaint()
{
	m_GLContext->makeCurrent();
	
	DetectOpenGLError();

	glDrawBuffer(GL_BACK);      // Set back buffer
	glClearColor((float)m_BkColor[0], (float)m_BkColor[1], (float)m_BkColor[2], (float)m_BkColor[3]);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); // Clear contents

}

void COpenGL::OnEndPaint()
{
	glFinish();
	m_GLContext->swapBuffers();
	m_GLContext->doneCurrent();
}



void COpenGL::OnDestroy()
{
	
}

/* ================== OpenGL Basic Operations ================== */
void COpenGL::Project(Coord pos, Coord& eye)
{
	m_GLContext->makeCurrent();

	GLint viewport[4];
	glGetIntegerv(GL_VIEWPORT, viewport); //Retrieves the viewport

	GLdouble modelview[16];
	glGetDoublev(GL_MODELVIEW_MATRIX, modelview); //Retrieve the matrix

	GLdouble projection[16];
	glGetDoublev(GL_PROJECTION_MATRIX, projection);
	gluProject(pos[0], pos[1], pos[2], modelview, projection, viewport, &(eye[0]), &(eye[1]), &(eye[2]));

	m_GLContext->doneCurrent();

}
void COpenGL::UnProject(Coord eye, Coord& pos)
{
	m_GLContext->makeCurrent();

	GLint viewport[4];
	glGetIntegerv(GL_VIEWPORT, viewport); //Retrieves the viewport

	GLdouble modelview[16];
	glGetDoublev(GL_MODELVIEW_MATRIX, modelview); //Retrieve the matrix

	GLdouble projection[16];
	glGetDoublev(GL_PROJECTION_MATRIX, projection);
	gluUnProject(eye[0], eye[1], eye[2], modelview, projection, viewport, &(pos[0]), &(pos[1]), &(pos[2]));

	m_GLContext->doneCurrent();
}
void COpenGL::GetDepth(int x, int y, float& depth)
{
	m_GLContext->makeCurrent();

	GLint viewport[4]; 
	glGetIntegerv(GL_VIEWPORT, viewport); //Retrieves the viewport

	glReadPixels(x, viewport[3]-y, 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, &depth);

	m_GLContext->doneCurrent();
}
void COpenGL::GetWorldCoord(int x, int y, Coord& pos)
{
	m_GLContext->makeCurrent();

	GLint viewport[4]; 
	glGetIntegerv(GL_VIEWPORT, viewport); //Retrieves the viewport

	GLdouble modelview[16];
	glGetDoublev(GL_MODELVIEW_MATRIX, modelview); //Retrieve the matrix

	GLdouble projection[16];
	glGetDoublev(GL_PROJECTION_MATRIX, projection);

	float winX, winY, winZ; //Window Coordinate	
	winX = (float)x; //Store the x coord
	winY = (float)(viewport[3] - y); //From Windows Coord System to OpenGL Coord System

	glReadPixels((GLint) winX, (GLint) winY, 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, &winZ);
	int  ret_val = gluUnProject(winX, winY, winZ, modelview, projection, viewport, &(pos[0]), &(pos[1]), &(pos[2]));
	assert(ret_val == GL_TRUE);

// 	std::cout<<"\nGetWorldCoord " << x <<" " <<y << std::endl;
// 
// 	std::cout << "Viewport : \n" << viewport[0] <<" " << viewport[1] <<" " << viewport[2] << " " << viewport[3] << std::endl;
// 	std::cout << "Modelview: \n" ;
// 	for(size_t k=0; k<4; ++k)
// 	{
// 		for(size_t i=0; i<4; ++i) std::cout << modelview[k*4+i] << " " ;
// 		std::cout << std::endl;
// 
// 	}
// 	std::cout << "Projection: \n";
// 	for(size_t k=0; k<4; ++k)
// 	{
// 		for(size_t i=0; i<4; ++i) std::cout << projection[k*4+i] << " ";
// 		std::cout << std::endl;
// 	}
// 	std::cout << winX <<" " << winY << " " << winZ << std::endl;
// 	std::cout << "\n GetWorldCoord end \n";


	m_GLContext->doneCurrent();
}
void COpenGL::GetViewPort(int viewport[4])
{
	m_GLContext->makeCurrent();

	glGetIntegerv(GL_VIEWPORT, viewport); //Retrieves the viewport

	m_GLContext->doneCurrent();
}
CMatrix COpenGL::GetModelViewMatrix()
{
	return m_ModelViewMatrix;
}
void COpenGL::SetModelViewMatrix(CMatrix& mvmatrix)
{
	m_ModelViewMatrix = mvmatrix;
}
void COpenGL::SetModelView()
{
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	double TempMatrix[16];
	m_ModelViewMatrix.Matrix2Array(TempMatrix);
	glMultMatrixd(TempMatrix);
}
bool COpenGL::LoadGLMatrix(std::string filename)
{
	std::ifstream file(filename.c_str());

	if (!file)
	{
		return false;
	}

	GLdouble modelview[16];
	GLdouble projection[16];

	for (size_t i = 0; i < 16; i++)
	{
		file >> modelview[i];
	}
	for (size_t i = 0; i < 16; i++)
	{
		file >> projection[i];
	}

	file.close();

	m_ModelViewMatrix.Array2Matrix(modelview);

	m_GLContext->makeCurrent();

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glMultMatrixd(modelview);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glMultMatrixd(projection);

	m_GLContext->doneCurrent();

	return true;
}
bool COpenGL::SaveGLMatrix(std::string filename)
{
	m_GLContext->makeCurrent();

	GLdouble modelview[16];
	glGetDoublev(GL_MODELVIEW_MATRIX, modelview); //Retrieve the matrix

	GLdouble projection[16];
	glGetDoublev(GL_PROJECTION_MATRIX, projection);

	m_GLContext->doneCurrent();

	std::ofstream file(filename.c_str());
	if (!file)
	{
		return false;
	}

	for (size_t i = 0; i < 4; i++)
	{
		for (size_t j = 0; j < 4; j++)
		{
			file << modelview[i * 4 + j] << ' ';
		}
		file << std::endl;
	}

	for (size_t i = 0; i < 4; i++)
	{
		for (size_t j = 0; j < 4; j++)
		{
			file << projection[i * 4 + j] << ' ';
		}
		file << std::endl;
	}

	file.close();
	return true;
}

void COpenGL::Create2DTexture(int type)
{
	// set check texture here.
	this->OnBeginPaint();

	switch(type)
	{
	case 1:
		MakeCheckImage();
		break;
	case 2:
		MakeLineImage();
		break;
	case 3:
		MakeBoundaryImage();
		break;
	default:
		break;
	}

	//
	glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
	glGenTextures(1, &texName);

	glBindTexture(GL_TEXTURE_2D, texName);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, CHECK_IMAGE_WIDTH, CHECK_IMAGE_HEIGHT, 
		0, GL_RGBA, GL_UNSIGNED_BYTE, checkImage);

	this->OnEndPaint();
}
void COpenGL::MakeCheckImage()
{
	size_t i, j, c;

	for (i = 0; i < CHECK_IMAGE_WIDTH; i++)
	{
		for (j = 0; j < CHECK_IMAGE_HEIGHT; j++)
		{	
			c = (((((i + 4) & 0x8) == 0) ^ ((((j+4) & 0x8)) == 0))) * 255;

			checkImage[i][j][0] = (GLubyte) c;
			checkImage[i][j][1] = (GLubyte) c;
			checkImage[i][j][2] = (GLubyte) c;
			checkImage[i][j][3] = (GLubyte) 255;
		}
	}

	for (i = 0; i < CHECK_IMAGE_WIDTH; i++)
	{
		checkImage[i][0][0] = (GLubyte) 255;
		checkImage[i][0][1] = (GLubyte) 0;
		checkImage[i][0][2] = (GLubyte) 0;
		checkImage[i][0][3] = (GLubyte) 255;

		checkImage[i][CHECK_IMAGE_HEIGHT - 1][0] = (GLubyte) 255;
		checkImage[i][CHECK_IMAGE_HEIGHT - 1][1] = (GLubyte) 0;
		checkImage[i][CHECK_IMAGE_HEIGHT - 1][2] = (GLubyte) 0;
		checkImage[i][CHECK_IMAGE_HEIGHT - 1][3] = (GLubyte) 255;
	}

	for (j = 0; j < CHECK_IMAGE_HEIGHT; j++)
	{
		checkImage[0][j][0] = (GLubyte) 255;
		checkImage[0][j][1] = (GLubyte) 0;
		checkImage[0][j][2] = (GLubyte) 0;
		checkImage[0][j][3] = (GLubyte) 255;

		checkImage[CHECK_IMAGE_WIDTH - 1][j][0] = (GLubyte) 255;
		checkImage[CHECK_IMAGE_WIDTH - 1][j][1] = (GLubyte) 0;
		checkImage[CHECK_IMAGE_WIDTH - 1][j][2] = (GLubyte) 0;
		checkImage[CHECK_IMAGE_WIDTH - 1][j][3] = (GLubyte) 255;
	}
}
void COpenGL::MakeLineImage()
{
	size_t i, j;

	// make a white background.
	for (i = 0; i < CHECK_IMAGE_WIDTH; i++)
	{
		for (j = 0; j < CHECK_IMAGE_HEIGHT; j++)
		{	
			checkImage[i][j][0] = (GLubyte) 255;
			checkImage[i][j][1] = (GLubyte) 255;
			checkImage[i][j][2] = (GLubyte) 255;
			checkImage[i][j][3] = (GLubyte) 255;
		}
	}

	// mask the lines.
	int distance = 4;
	for (i = distance - 1; i < CHECK_IMAGE_WIDTH; i+=distance)
	{
		for (j = 0; j < CHECK_IMAGE_HEIGHT; j++)
		{	
			checkImage[i][j][0] = (GLubyte) 255;
			checkImage[i][j][1] = (GLubyte) 0;
			checkImage[i][j][2] = (GLubyte) 0;
			checkImage[i][j][3] = (GLubyte) 255;
		}
	}
	for (j = distance - 1; j < CHECK_IMAGE_HEIGHT; j+=distance)
	{	
		for (i = 0; i < CHECK_IMAGE_WIDTH; i++)
		{
			checkImage[i][j][0] = (GLubyte) 0;
			checkImage[i][j][1] = (GLubyte) 0;
			checkImage[i][j][2] = (GLubyte) 255;
			checkImage[i][j][3] = (GLubyte) 255;
		}
	}
}
void COpenGL::MakeBoundaryImage()
{
	size_t i, j;

	// make a white background.
	for (i = 0; i < CHECK_IMAGE_WIDTH; i++)
	{
		for (j = 0; j < CHECK_IMAGE_HEIGHT; j++)
		{	
			checkImage[i][j][0] = (GLubyte) 255;
			checkImage[i][j][1] = (GLubyte) 255;
			checkImage[i][j][2] = (GLubyte) 255;
			checkImage[i][j][3] = (GLubyte) 255;
		}
	}

	// set boundary here.
	for (i = 0; i < CHECK_IMAGE_WIDTH; i++)
	{
		checkImage[i][0][0] = (GLubyte) 0;
		checkImage[i][0][1] = (GLubyte) 0;
		checkImage[i][0][2] = (GLubyte) 0;
		checkImage[i][0][3] = (GLubyte) 0;

		checkImage[i][CHECK_IMAGE_HEIGHT - 1][0] = (GLubyte) 0;
		checkImage[i][CHECK_IMAGE_HEIGHT - 1][1] = (GLubyte) 0;
		checkImage[i][CHECK_IMAGE_HEIGHT - 1][2] = (GLubyte) 0;
		checkImage[i][CHECK_IMAGE_HEIGHT - 1][3] = (GLubyte) 0;
	}

	for (j = 0; j < CHECK_IMAGE_HEIGHT; j++)
	{
		checkImage[0][j][0] = (GLubyte) 0;
		checkImage[0][j][1] = (GLubyte) 0;
		checkImage[0][j][2] = (GLubyte) 0;
		checkImage[0][j][3] = (GLubyte) 0;

		checkImage[CHECK_IMAGE_WIDTH - 1][j][0] = (GLubyte) 0;
		checkImage[CHECK_IMAGE_WIDTH - 1][j][1] = (GLubyte) 0;
		checkImage[CHECK_IMAGE_WIDTH - 1][j][2] = (GLubyte) 0;
		checkImage[CHECK_IMAGE_WIDTH - 1][j][3] = (GLubyte) 0;
	}
}
