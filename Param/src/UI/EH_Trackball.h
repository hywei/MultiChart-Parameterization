#ifndef EVENT_HANDLER_TRACKBALL_H
#define EVENT_HANDLER_TRACKBALL_H

#include "EventHandlerInterface.h"
#include "../Numerical/matrix.h"
#include "../Numerical/Trackball.h"
#include <set>

class EH_Trackball : public EventHandlerInterface
{
public:
	EH_Trackball();
	~EH_Trackball();

public:
	void Init(QGLViewer* _qtViewer);
	bool OnLButtonDown(unsigned int nFlags, int point_x, int point_y);
	bool OnMouseMove(unsigned int nFlags, int point_x, int point_y);
	bool OnMouseWheel(unsigned int nFlags, short zDelta, int point_x, int point_y);
	handler_key HandlerKey() {return "TrackballHandler";} 
	bool IsMouseNeedHandle(MouseMode mouse_mode);


private:
	 CMatrix      m_ViewMtx;     // Model/View matrix
	 Trackball    m_Trackball;   // Track ball structure
	 std::pair<int, int> m_OldPoint;
	 std::set<MouseMode> m_ValidMouseMode;  
};
#endif
