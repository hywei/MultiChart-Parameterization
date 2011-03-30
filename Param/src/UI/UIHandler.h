#ifndef UI_HANDLER_H
#define UI_HANDLER_H

#include "EventHandlerInterface.h"

#include <vector>
using namespace std;

class CUIHandler
{
public:
	CUIHandler();
	~CUIHandler();

public:
	void Init(QGLViewer* _glViewer);
	bool OnLButtonDown(unsigned int nFlags, int point_x, int point_y);
	bool OnLButtonUp(unsigned int nFlags, int point_x, int point_y);
	bool OnRButtonDown(unsigned int nFlags, int point_x, int point_y);
	bool OnRButtonUp(unsigned int nFlags, int point_x, int point_y);
	bool OnMouseMove(unsigned int nFlags, int point_x, int point_y);
	bool OnMouseWheel(unsigned int nFlags, short zDelta, int point_x, int point_y);
	bool OnKeyDown(unsigned int nChar);

public:
	MouseMode& GetMouseMode() {return m_MouseMode;}

private:
	MouseMode m_MouseMode;
	KeyMode   m_KeyMode;
	vector<EventHandlerInterface*>  m_EH_Layer;
};

#endif