
#include "UIHandler.h"
#include "EH_Trackball.h"
#include "EH_Param.h"
#include "EH_Cmd.h"
#include <cassert>


//////////////////////////////////////////////////////////////////////
// CUIHandler Construction/Destruction
//////////////////////////////////////////////////////////////////////
CUIHandler::CUIHandler()
{
	m_EH_Layer.push_back(new EH_Trackball);
	m_EH_Layer.push_back(new EH_Param);
}
CUIHandler::~CUIHandler()
{
	for (size_t i = 0; i < m_EH_Layer.size(); i++)
	{
		delete m_EH_Layer[i];
	}
}

//////////////////////////////////////////////////////////////////////
// CUIHandler event handler methods
//////////////////////////////////////////////////////////////////////
void CUIHandler::Init(QGLViewer* _glViewer)
{
	assert(_glViewer != NULL);

	for (size_t i = 0; i < m_EH_Layer.size(); i++)
	{
		m_EH_Layer[i]->Init(_glViewer);
	}
}
bool CUIHandler::OnLButtonDown(unsigned int nFlags, int point_x, int point_y)
{
	bool ret_val = true;
	for (size_t i = 0; i < m_EH_Layer.size(); i++)
	{
		EventHandlerInterface* eh = m_EH_Layer[i];
		
		if (eh->IsMouseNeedHandle(m_MouseMode))
		{
			ret_val &= eh->OnLButtonDown(nFlags, point_x, point_y);
		}
	}

	return ret_val;
}
bool CUIHandler::OnLButtonUp(unsigned int nFlags, int point_x, int point_y)
{
	bool ret_val = true;
	for (size_t i = 0; i < m_EH_Layer.size(); i++)
	{
		EventHandlerInterface* eh = m_EH_Layer[i];

		if (eh->IsMouseNeedHandle(m_MouseMode))
		{
			ret_val &= eh->OnLButtonUp(nFlags, point_x, point_y);
		}
	}

	return ret_val;
}
bool CUIHandler::OnRButtonDown(unsigned int nFlags, int point_x, int point_y)
{
	bool ret_val = true;
	for (size_t i = 0; i < m_EH_Layer.size(); i++)
	{
		EventHandlerInterface* eh = m_EH_Layer[i];

		if (eh->IsMouseNeedHandle(m_MouseMode))
		{
			ret_val &= eh->OnRButtonDown(nFlags, point_x, point_y);
		}
	}

	return ret_val;
}
bool CUIHandler::OnRButtonUp(unsigned int nFlags, int point_x, int point_y)
{
	bool ret_val = true;
	for (size_t i = 0; i < m_EH_Layer.size(); i++)
	{
		EventHandlerInterface* eh = m_EH_Layer[i];

		if (eh->IsMouseNeedHandle(m_MouseMode))
		{
			ret_val &= eh->OnRButtonUp(nFlags, point_x, point_y);
		}
	}

	return ret_val;
}
bool CUIHandler::OnMouseMove(unsigned int nFlags, int point_x, int point_y)
{
	bool ret_val = true;
	for (size_t i = 0; i < m_EH_Layer.size(); i++)
	{
		EventHandlerInterface* eh = m_EH_Layer[i];

		if (eh->IsMouseNeedHandle(m_MouseMode))
		{
			ret_val &= eh->OnMouseMove(nFlags, point_x, point_y);
		}
	}

	return ret_val;
}
bool CUIHandler::OnMouseWheel(unsigned int nFlags, short zDelta, int point_x, int point_y)
{
	bool ret_val = true;
	for (size_t i = 0; i < m_EH_Layer.size(); i++)
	{
		EventHandlerInterface* eh = m_EH_Layer[i];

		ret_val &= eh->OnMouseWheel(nFlags, zDelta, point_x, point_y);
	}

	return ret_val;
}
bool CUIHandler::OnKeyDown(unsigned int nChar)
{
	bool ret_val = true;
	for (size_t i = 0; i < m_EH_Layer.size(); i++)
	{
		EventHandlerInterface* eh = m_EH_Layer[i];

		if (eh->IsKeyNeedHandle(m_KeyMode))
		{
			ret_val &= eh->OnKeyDown(nChar);
		}
	}

	return ret_val;
}