#ifndef OPENMESHAPPS_QGLViewer_HH
#define OPENMESHAPPS_QGLViewer_HH


//== INCLUDES =================================================================

#include "../Common/BasicDataType.h"
#include <QtOpenGL>
#include <boost/shared_ptr.hpp>

//== FORWARD DECLARATIONS =====================================================

class QMenu;
class MeshModel;
class COpenGL;
class CUIHandler;
class Utility;

namespace PARAM
{
    class Parameter;
	class ParamDrawer;
}

//== CLASS DEFINITION =========================================================


class QGLViewer : public QGLWidget
{

	Q_OBJECT

private:
	// 

public:
	typedef QGLWidget Super;

	// Default constructor.
	QGLViewer( QWidget* _parent=0 );
	QGLViewer( QGLFormat& _fmt, QWidget* _parent=0 );
	// Destructor.
	virtual ~QGLViewer();


public:
	void getBoundingSphere(Coord& center, double& radius);

	void EmitSelectVertexSignal()
	{
		emit select_vertex();
	}

signals:
	void select_vertex(void);

public slots:
	void loadMeshModel();
	void saveMeshModel();
	int loadTextureImage();
	int loadQuadFile();
	
	void mouseSpin();
	void mouseMove();
	void mouseZoom();

	void CreateSquareTexture();
	void CreateLineTexture();
	void CreateBoundaryTexture();
    
	void SolveParameter();

	void SaveFaceTexCoord();
	void LoadFaceTexCoord();

	void SetParamDrawerSelectPatchMode();
	void SetParamDrawerSelectVertMode();
	void SetParamDrawerCorrespondMode();

    void SetChartInitValue(const std::vector<double>&);
    void ChartOptimization();
private slots:
	void RenderingSolidSmooth();
	void RenderingSolidFlat();
	void RenderingSolidWireframe();
	void RenderingTransparent();
	void RenderingWireframe();
	void RenderingVertices();
	void RenderingVertexColor();
	void RenderingFaceColor();
	void RenderingTexture();
	void RenderingBoundingBox();
	void RenderingParamTexture();


private: // inherited

	// initialize OpenGL states (triggered by Qt)
	void initializeGL();

	// draw the scene (triggered by Qt)
	void paintGL();

	// handle resize events (triggered by Qt)
	void resizeGL( int w, int h );

private:
	void init();
	void setDefaultMaterial();
	void setDefaultLight();

	void createPopMenu();

	int CreateTexture(const std::string& texture_file_name);


protected:

	// Qt mouse events
	virtual void mousePressEvent( QMouseEvent* );
	virtual void mouseReleaseEvent( QMouseEvent* );
	virtual void mouseMoveEvent( QMouseEvent* );
	virtual void wheelEvent( QWheelEvent* );
	virtual void keyPressEvent( QKeyEvent* );

private:
	// popup menu for draw mode selection
	QMenu*	popup_menu_;

	GLuint texture[2];
public:
	//
	boost::shared_ptr<MeshModel> p_mesh;
	boost::shared_ptr<COpenGL> p_opengl;
	boost::shared_ptr<CUIHandler> p_UIHander;
	boost::shared_ptr<Utility> p_util;
	boost::shared_ptr<PARAM::Parameter> p_param;
	boost::shared_ptr<PARAM::ParamDrawer> p_param_drawer;
};


//=============================================================================
#endif // OPENMESHAPPS_QGLViewer_HH
//=============================================================================
