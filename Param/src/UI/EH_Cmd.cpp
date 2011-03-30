#include "EH_Cmd.h"
#include "../Common/BasicDataType.h"
#include "../MainWindow/QGLViewer.h"

//////////////////////////////////////////////////////////////////////
// EH_Cmd Construction/Destruction
//////////////////////////////////////////////////////////////////////
EH_Cmd::EH_Cmd()
{
	this->m_MouseMode = MOUSE_MODE_UNDEFINED;
	this->m_KeyMode = KEY_MODE_UNDEFINED;
}
EH_Cmd::~EH_Cmd()
{
	this->glViewer = NULL;
}
//////////////////////////////////////////////////////////////////////
// EH_Cmd event handler methods
//////////////////////////////////////////////////////////////////////
void EH_Cmd::Init(QGLViewer* _glViewer)
{
	this->glViewer = _glViewer; 
}
bool EH_Cmd::OnKeyDown(unsigned int nChar)
{
	return true;
}
bool EH_Cmd::IsKeyNeedHandle(KeyMode key_mode)
{
	this->m_KeyMode = key_mode;
	return true;
}