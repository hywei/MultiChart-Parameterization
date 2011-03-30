/*_____________________________________________________
 |
 |                                          
 |  文件：Trackball.cpp
 |
*/
	
#include <math.h>		
#include <assert.h>
#include "Trackball.h"


Trackball::Trackball(void)
{
	m_iWidth=m_iHeight=400;	
	m_iCenterX=m_iCenterY=200;
	m_iRadius=200;
	m_iMouseX=m_iMouseY=1;

	// 零矢量
	m_start[0]=0;
	m_start[1]=0;
	m_start[2]=0;
}

Trackball::~Trackball(void)
{

}

// 跟踪球的中心为窗口的中心
void Trackball::SetTrackWindow(int iWidth, int iHeight)
{
	m_iWidth   = iWidth;
	m_iHeight  = iHeight;
	m_iCenterX = m_iWidth/2;
	m_iCenterY = m_iHeight/2;
	m_iRadius  = (m_iWidth > m_iHeight) ? m_iCenterY : m_iCenterX ;
}

// 跟踪球的中心为窗口坐标(centerX, centerY) 
void Trackball::SetTrackWindow(int iWidth, int iHeight, int iCenterX, int iCenterY)
{
	m_iWidth   = iWidth;
	m_iHeight  = iHeight;
	m_iCenterX = iCenterX;
	m_iCenterY = iCenterY;
	m_iRadius  = (m_iWidth > m_iHeight) ? m_iCenterY : m_iCenterX ;
}

// 将鼠标(mouseX, mouseY)投影到球面一点 P
// 函数值即为球心到P的矢量
Coord Trackball::MouseToCoord(int iMouseX, int iMouseY)
{
	Coord V(0,0,0);

	if(m_iCenterX >= m_iWidth/2)  
		V[0]=float(iMouseX-m_iCenterX)/float(m_iCenterX);
	else
		V[0]=float(iMouseX-m_iCenterX)/float(m_iWidth-m_iCenterX);

	if(m_iCenterY >= m_iHeight/2)  
		V[1]=float(m_iCenterY-iMouseY)/float(m_iCenterY);
	else
		V[1]=float(m_iCenterY-iMouseY)/float(m_iHeight-m_iCenterY);

	float d=float(sqrt(V[0] * V[0] + V[1] * V[1]));
	V[2]=float( cos((3.14159265 / 2.0) * ((d < 1.0) ? d : 1.0)) );
	V.normalize();
	return V;
}

// 输入刚按下鼠标时的坐标 (mouseX, mouseY)
void Trackball::Start(int iMouseX, int iMouseY)
{
	m_iMouseX=iMouseX;
	m_iMouseY=iMouseY;
	m_start=MouseToCoord(m_iMouseX, m_iMouseY);
}

// 鼠标移动 (dx,dy), 计算旋转轴 axis 与角度 angle 
void Trackball::Tracking(int iDx, int iDy, Coord *axis, float *fAngle)
{
	m_iMouseX+=iDx;
	m_iMouseY+=iDy;
	m_stop=MouseToCoord(m_iMouseX, m_iMouseY);
	*axis=cross(m_start, m_stop); // cross product

	const float ANGLE_SCALE=2.0f;	// 任意给定一个角度放大系数
	*fAngle=float(ANGLE_SCALE*acos(dot(m_start, m_stop))/3.14159*180.0);
	
	// 控制输出的最小角度，否则 Rotation 将出错误
	const float	DELTA_ANGLE=0.00001f;
	if(*fAngle < DELTA_ANGLE)	*fAngle=DELTA_ANGLE;

	m_start=m_stop; // 不要忘记
}

/*-------------------------------------------------------
	复杂的跟踪球算法可以参考 
	蔡文立，《科学计算可视化算法与系统》 P232 - P233
  -------------------------------------------------------*/