#include "EH_Trackball.h"
#include "../MainWindow/QGLViewer.h"
#include "../Common/BasicDataType.h"
#include "../OpenGL/GLElement.h"

#include <iostream>

#define MOUSE_ZOOM_RATIO_LARGE      0.1
#define MOUSE_ZOOM_RATIO_MIDDLE     0.01
#define MOUSE_ZOOM_RATIO_SMALL      0.001
 
//////////////////////////////////////////////////////////////////////
// EH_Trackball Construction/Destruction
//////////////////////////////////////////////////////////////////////
EH_Trackball::EH_Trackball()
{
	this->m_MouseMode = MOUSE_MODE_UNDEFINED;
	this->m_KeyMode = KEY_MODE_UNDEFINED;

	//
	m_ValidMouseMode.insert(MOUSE_MODE_MOVE);
	m_ValidMouseMode.insert(MOUSE_MODE_SPIN);
	m_ValidMouseMode.insert(MOUSE_MODE_ZOOM);
}
EH_Trackball::~EH_Trackball()
{
	this->glViewer = NULL;
}

//////////////////////////////////////////////////////////////////////
// EH_Trackball event handler methods
//////////////////////////////////////////////////////////////////////
void EH_Trackball::Init(QGLViewer* _glViewer)
{
	this->glViewer = _glViewer;
}
bool EH_Trackball::OnLButtonDown(unsigned int nFlags, int point_x, int point_y)
{
	switch(m_MouseMode) 
	{
	case MOUSE_MODE_MOVE:
		break;

	case MOUSE_MODE_SPIN:
		{
			GLint viewport[4]; 
			glViewer->p_opengl->GetViewPort(viewport);
			
			m_Trackball.SetTrackWindow(viewport[2], viewport[3]);
			m_Trackball.Start(point_x, point_y);
		}
		break;

	case MOUSE_MODE_ZOOM:
		break;
	}

	m_OldPoint = std::make_pair(point_x, point_y);
	return true;
}
bool EH_Trackball::OnMouseMove(unsigned int nFlags, int point_x, int point_y)
{
	// LButton is pressed (dragging)
	if((nFlags&Qt::LeftButton) != Qt::LeftButton)   
	{
		return true;
	}

	Coord Center;
	double radius = 0;
	glViewer->getBoundingSphere(Center, radius);
	
	m_ViewMtx = glViewer->p_opengl->GetModelViewMatrix();
	

	switch(m_MouseMode) {
	case MOUSE_MODE_MOVE:
		{
			double dx,dy,dz;
			// Method thinking in Eye Coordinate System
			Coord CenterAfter, NewCenter;
			glViewer->p_opengl->Project(Center, CenterAfter);

			CenterAfter[0] += point_x - m_OldPoint.first;
			CenterAfter[1] += m_OldPoint.second - point_y;
			glViewer->p_opengl->UnProject(CenterAfter, NewCenter);
		
			Coord CenterView = Center;
			Coord NewCenterView = NewCenter;
			m_ViewMtx.TransformVector(CenterView);
			m_ViewMtx.TransformVector(NewCenterView);

			dx = NewCenterView[0] - CenterView[0];
			dy = NewCenterView[1] - CenterView[1];
			dz = NewCenterView[2] - CenterView[2];
			// 			fprintf(stdout, "Translate(%f, %f, %f)\n", dx, dy, dz);

			CMatrix TempMatrix(4);
			TempMatrix.SetTrans2Matrix(Coord(dx, dy, dz));
			m_ViewMtx = TempMatrix*m_ViewMtx;			

		}
		break;

	case MOUSE_MODE_SPIN:
		{		
			Coord pAxis;
			float pAngle=0;
			m_Trackball.Tracking((point_x-m_OldPoint.first), (point_y-m_OldPoint.second), &pAxis, &pAngle);

			CMatrix PreTransMtx(4), PostTransMtx(4), RotMtx(4);

			Coord CenterAfter = Center;
			m_ViewMtx.TransformVector(CenterAfter);

			PreTransMtx.SetTrans2Matrix(-CenterAfter);
			PostTransMtx.SetTrans2Matrix(CenterAfter);

			// Constrain Rotation Axis
			if((nFlags&Qt::CTRL) && (nFlags&Qt::SHIFT))	// Around Z Axis
			{
				if(pAxis[0] > 0.0)
					pAxis.setCoords(0.0,0.0,1.0);
				else
					pAxis.setCoords(0.0,0.0,-1.0);
			}
			else 
			{
				if(nFlags&Qt::CTRL)	// Around X Axis
				{
					if(pAxis[0] > 0.0)
						pAxis.setCoords(1.0,0.0,0.0);
					else
						pAxis.setCoords(-1.0,0.0,0.0);
				}
				else if(nFlags&Qt::SHIFT)// Around Y Axis
				{
					if(pAxis[0] > 0.0)
						pAxis.setCoords(0.0,1.0,0.0);
					else
						pAxis.setCoords(0.0,-1.0,0.0);
				}
			}

			RotMtx.SetRotation2Matrix(pAngle, pAxis);
			m_ViewMtx = PostTransMtx*RotMtx*PreTransMtx*m_ViewMtx;
		}        
		break; 
	}

	m_OldPoint = std::make_pair(point_x, point_y);
	glViewer->p_opengl->SetModelViewMatrix(m_ViewMtx);

	return true;
}
bool EH_Trackball::OnMouseWheel(unsigned int nFlags, short zDelta, int point_x, int point_y)
{

	double ratio;
	if(nFlags&Qt::CTRL)
		ratio = MOUSE_ZOOM_RATIO_LARGE;
	else if(nFlags&Qt::SHIFT)
		ratio = MOUSE_ZOOM_RATIO_SMALL;
	else
		ratio = MOUSE_ZOOM_RATIO_MIDDLE;

	double scale = (zDelta > 0) ? (1.0 + ratio) : (1.0 - ratio);
	glViewer->p_opengl->m_Projection.Zoom(scale);

	glViewer->p_opengl->OnResize(glViewer->width(), glViewer->height());

	return true;
}
bool EH_Trackball::IsMouseNeedHandle(MouseMode mouse_mode)
{
	this->m_MouseMode = mouse_mode;
	return m_ValidMouseMode.find(mouse_mode) != m_ValidMouseMode.end();
}