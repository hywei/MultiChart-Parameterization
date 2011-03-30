#include "GLContext.h"

GLContext::GLContext(const QGLFormat &fmt, QPaintDevice *pDevice):
		QGLContext(fmt, pDevice)
{
	m_nEnterCounter=0;
}

GLContext::~GLContext() {}

bool GLContext::OnCreate()
{
	return create();
}

void GLContext::OnDestroy()
{
}

void GLContext::glEnter()
{
	makeCurrent();
}

void GLContext::glLeave()
{
	m_nEnterCounter--;
	if(m_nEnterCounter==0)
		doneCurrent();
}

void GLContext::SwapBuffers()
{
	swapBuffers();
}


