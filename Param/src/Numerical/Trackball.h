/*_____________________________________________________
 |
 |  文件：Trackball.h
 |
 |_____________________________________________________*/




#if !defined(TRACKBALL_H)
#define		 TRACKBALL_H

#include "../Common/BasicDataType.h"

class Trackball
{
public:
	Trackball(void);
   ~Trackball(void);

public:
	// 跟踪球的中心为窗口的中心
	void SetTrackWindow(int iWidth, int iHeight);

	// 跟踪球的中心为窗口坐标(centerX, centerY) 
	void SetTrackWindow(int iWidth, int iHeight, int iCenterX, int iCenterY);

	// 输入刚按下鼠标时的坐标 (mouseX, mouseY)
	void Start(int iMouseX, int iMouseY);

	// 鼠标移动 (dx,dy), 计算旋转轴 axis 与角度 angle 
	void Tracking(int iDx, int iDy, Coord *axis, float *fAngle);

protected:

	// 将鼠标(mouseX, mouseY)投影到球面一点 P
	// 函数值即为球心到P的矢量
	Coord MouseToCoord(int iMouseX, int iMouseY);

	int    m_iWidth,   m_iHeight;	// 窗口的宽高
	int    m_iCenterX, m_iCenterY;	// 跟踪球中心
	int    m_iRadius;				// 跟踪球半径
	int    m_iMouseX,  m_iMouseY;	// 鼠标位置

	Coord m_start;	// 刚按下时的鼠标投影矢量
	Coord m_stop;	// 移动(dx,dy)时的鼠标投影矢量
};

#endif // TRACKBALL_H
