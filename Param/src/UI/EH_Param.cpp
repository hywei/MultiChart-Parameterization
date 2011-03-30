#include "EH_Param.h"
#include "../MainWindow/QGLViewer.h"
#include "../OpenGL/GLElement.h"
#include "../Param/ParamDrawer.h"
#include <iostream>

EH_Param::EH_Param()
{
	this->m_MouseMode = MOUSE_MODE_UNDEFINED;
	this->m_KeyMode = KEY_MODE_UNDEFINED;

	m_ValidMouseMode.insert(MOUSE_MODE_SELECT_VERTEX);
	m_ValidMouseMode.insert(MOUSE_MODE_SELECT_PATCH);
}

EH_Param::~EH_Param()
{
	glViewer = NULL;
}

void EH_Param::Init(QGLViewer* _glViewer)
{
	glViewer = _glViewer;
}

bool EH_Param::OnLButtonDown(unsigned int nFlags, int point_x, int point_y)
{
//	std::cout << "LButtonDown : " << nFlags << " " << point_x << " " << point_y << std::endl;
	Coord hit_coord;
	switch(m_MouseMode)
	{
	case MOUSE_MODE_SELECT_VERTEX:
		glViewer->p_opengl->GetWorldCoord(point_x, point_y, hit_coord);
		glViewer->p_param_drawer->SetSelectedVertCoord(hit_coord);
//         std::cout << "Left Button : " << hit_coord[0] << ' ' << hit_coord[1] <<
//             ' ' << hit_coord[2] << std::endl;
		glViewer->EmitSelectVertexSignal();
		break;

	case MOUSE_MODE_SELECT_PATCH:
		glViewer->p_opengl->GetWorldCoord(point_x, point_y, hit_coord);
		glViewer->p_param_drawer->SetSelectedPatchCoord(hit_coord);
		break;


	default:
		break;
	}
	return true;
}

bool EH_Param::OnMouseMove(unsigned int nFlags, int point_x, int point_y)
{
//	std::cout <<"Mouse Move : " << nFlags << " " << point_x << " " << point_y << std::endl;
	Coord hit_coord;
	switch(m_MouseMode)
	{
	case MOUSE_MODE_SELECT_VERTEX:
		glViewer->p_opengl->GetWorldCoord(point_x, point_y, hit_coord);
		glViewer->p_param_drawer->SetSelectedVertCoord(hit_coord);
//         std::cout << "Mouse Move : " << hit_coord[0] << ' ' << hit_coord[1] <<
//             ' ' << hit_coord[2] << std::endl;
		glViewer->EmitSelectVertexSignal();
		break;

	default:
		break;
	}
	return true;
}

bool EH_Param::IsMouseNeedHandle(MouseMode mouse_mode)
{
	this->m_MouseMode = mouse_mode;
	return m_ValidMouseMode.find(mouse_mode) != m_ValidMouseMode.end();
}
