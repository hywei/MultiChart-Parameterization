#ifndef EVENT_INTERFACE_H
#define EVENT_INTERFACE_H

#include <string>

typedef std::string handler_key;

class QGLViewer;

typedef enum { MOUSE_MODE_UNDEFINED,
	MOUSE_MODE_MOVE, MOUSE_MODE_SPIN, MOUSE_MODE_ZOOM,
	MOUSE_MODE_SELECT_VERTEX,
	MOUSE_MODE_SELECT_PATCH
} MouseMode;

typedef enum {
	KEY_MODE_UNDEFINED
} KeyMode;

class EventHandlerInterface
{
public:
	virtual ~EventHandlerInterface(){};

public:
	virtual void Init(QGLViewer* _qtViewer) = 0;
	virtual bool OnLButtonDown(unsigned int nFlags, int point_x, int point_y) {return true;}
	virtual bool OnLButtonUp(unsigned int nFlags, int point_x, int point_y) {return true;}
	virtual bool OnRButtonDown(unsigned int nFlags, int point_x, int point_y) {return true;}
	virtual bool OnRButtonUp(unsigned int nFlags, int point_x, int point_y) {return true;}
	virtual bool OnMouseMove(unsigned int nFlags, int point_x, int point_y) {return true;}
	virtual bool OnMouseWheel(unsigned int nFlags, short zDelta, int point_x, int point_y) {return true;}
	virtual bool OnKeyDown(unsigned int nChar) {return true;}
	virtual handler_key HandlerKey() {return "";}
	virtual bool IsMouseNeedHandle(MouseMode mouse_mode) {return false;}
	virtual bool IsKeyNeedHandle(KeyMode key_mode) {return false;}

protected:
	MouseMode m_MouseMode;
	KeyMode m_KeyMode;

	QGLViewer* glViewer;
};

#endif
