#ifndef GLCONTEXT_H
#define GLCONTEXT_H

#include <QGLContext>

class GLContext : public QGLContext
{
public:
	GLContext(const QGLFormat& fmt, QPaintDevice* pDevice);

	virtual ~GLContext();

public:
	bool OnCreate();
	void OnDestroy();

	void glEnter();
	void glLeave();

	void SwapBuffers();

private:
	int m_nEnterCounter;
};

#endif // GLCONTEXT_H
