#ifndef EH_PARAM_H_
#define EH_PARAM_H_

#include "EventHandlerInterface.h"
#include <set>

class EH_Param : public EventHandlerInterface
{
public:
	EH_Param();
	~EH_Param();

public:
	void Init(QGLViewer* _glViewer);
	bool OnLButtonDown(unsigned int nFlags, int point_x, int point_y);
	bool OnMouseMove(unsigned int nFlags, int point_x, int point_y);

	handler_key HandleKey() { return "ParameterHandler"; }
	bool IsMouseNeedHandle(MouseMode mouse_mode);

private:
	std::set<MouseMode> m_ValidMouseMode;
};

#endif