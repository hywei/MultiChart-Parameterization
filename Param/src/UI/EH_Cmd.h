#ifndef EVENT_HANDLER_CMD_H
#define EVENT_HANDLER_CMD_H

#include "EventHandlerInterface.h"

class EH_Cmd : public EventHandlerInterface
{
public:
	EH_Cmd();
	~EH_Cmd();

public:
	void Init(QGLViewer* _glViewer);
	handler_key HandlerKey() {return "CmdHandler";} 
	bool OnKeyDown(unsigned int nChar);
	bool IsKeyNeedHandle(KeyMode key_mode);
};

#endif