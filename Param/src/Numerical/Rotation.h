/*_____________________________________________________
 |
 |	Intra3D Layer1 :  Algebra 构件组成之一
 |                                          
 |  文件：Rotation.h
 |
 |  功能：四元组旋转运算
 |
 |  开发：林锐 ，1999/01/10
 |
 |	源程序测试：进行了单步跟踪测试
 |_____________________________________________________*/

/*_____________________________________________________
 |  
 |	基于 OpenGL 的交互式三维图形软件开发工具
 |
 |	**   Intra3D 2.0 for Windows 9X/NT  **
 |
 |	          著作权人：林锐 
 |	
 |	浙江大学 CAD&CG 国家重点实验室 (310027)
 |_____________________________________________________*/


#if !defined(ROTATION_H)
#define		 ROTATION_H

#include "../Common/BasicDataType.h"

/*______________________ ROTATION 通过 QUATERNION 实现 _____________________*/


// 定义旋转 ROTATION 变换, ROTATION 归一化是指 axis 归一化
// 用 ROTATION 结构表示旋转时，(axis, angle) 与 -(axis, angle) 等价

class ROTATION
{
public:
	Coord axis ;	// Unit Vector  
	double  angle;	// 0 - 360 degree
	// 构造归一化的 ROTATION 
	ROTATION(Coord Axis, double fAngle=0.0f)
	{ axis=Axis; angle=fAngle; axis.normalize(); }
	// 构造归一化的 ROTATION 
	ROTATION(double fX=0.0f, double fY=0.0f, double fZ=1.0f, double fAngle=0.0f)
	{ axis[0]=fX; axis[1]=fY; axis[2]=fZ; angle=fAngle; axis.normalize(); }

	
};

// ROTATION 相乘，返回值 = R1 * R2  (先执行 R1 旋转，后执行 R2 旋转)
// 要求 R1 R2 先归一化
ROTATION RotationMultiply( ROTATION R1,  ROTATION R2);

ROTATION operator *      ( ROTATION R1,  ROTATION R2);

// V 经过 R 变换后为 V', 函数返回 V'
// 要求 R 先归一化

Coord VectorTransform(const Coord V, const ROTATION R);

// 将 ROTATION 结构表示成 Matrix
// 要求 R 先归一化
void RotationToMatrix(double M[16], const ROTATION R);


/*_______________ 不想了解四元组细节的程序员不必往下看 _____________________*/



// 定义四元组 QUATERNION
// 用 QUATERNION 结构表示旋转时，(x,y,z,w) 与 -(x,y,z,w) 等价

class QUATERNION
{
public:
	double x, y, z, w;
	QUATERNION(double fX=0.0f, double fY=0.0f, double fZ=0.0f, double fW=1.0f)
	{ x=fX; y=fY; z=fZ; w=fW; }
};


// 四元组求模，返回值 = |A|
double QuaternionMagnitude(const QUATERNION A);

// 四元组归一化
// 如果 |A|=0，输出值 = A；输出值 = A/(|A|)
void QuaternionNormalize(QUATERNION *A);

// 四元组求逆
// 如果 |A|=0，输出值 = A；否则输出值 = A 的逆
void QuaternionInverse(QUATERNION *A);

// 四元组共扼，输出值 = A 的共扼
void QuaternionConjugate(QUATERNION *A);

// 四元组相加，返回值 = A + B
QUATERNION QuaternionAdd(const QUATERNION A, const QUATERNION B);

QUATERNION operator +   (const QUATERNION A, const QUATERNION B);

// 四元组相减，返回值 = A - B
QUATERNION QuaternionSub(const QUATERNION A, const QUATERNION B);

QUATERNION operator -   (const QUATERNION A, const QUATERNION B);

// 四元组缩放，返回值 = s * A 
QUATERNION QuaternionScale(const QUATERNION A, const double s);

QUATERNION operator *     (const QUATERNION A, const double s);

QUATERNION QuaternionScale(const double s, const QUATERNION A);

QUATERNION operator *     (const double s, const QUATERNION A);

// 四元组相乘，返回值 = A * B
QUATERNION QuaternionMultiply(const QUATERNION A, const QUATERNION B);

QUATERNION operator *        (const QUATERNION A, const QUATERNION B);  

// Spherical Linear Interpolation Between two Quaternions, 0 < t < 1
QUATERNION QuaternionSlerp(const QUATERNION from, const QUATERNION to, double t);

/*_____________________________________________________
 |                                                     
 |   为提高计算性能，以下变换函数均假定
 |
 |   输入参数 R Q V 已经进行了 Normalize 处理    
 |_____________________________________________________*/


/*------------------------------------------------------------------*/
/*----------------------  QUATERNION — ROTATION  ------------------*/

// 将 ROTATION 结构表示成 QUATERNION
QUATERNION RotationToQuaternion(const ROTATION R);

// 将 QUATERNION 结构表示成 ROTATION
ROTATION QuaternionToRotation(const QUATERNION Q);

/*------------------------------------------------------------------*/
/*----------------------  QUATERNION — Matrix  --------------------*/

// 将 QUATERNION 结构表示成 Matrix
void QuaternionToMatrix(double M[16], const QUATERNION Q);
void QuaternionToMatrix(double R[3][3], const QUATERNION quat);

// 将 Matrix 结构表示成 QUATERNION
// 要求 Matrix 是一种旋转矩阵,否则不能得到正确结果
void MatrixToQuaternion(QUATERNION *Q, const double M[16]);
void MatrixToQuaternion(QUATERNION *quat, const double R[3][3]);
void MatrixToQuaternion2(QUATERNION& quat, const double R[3][3]);
/*------------------------------------------------------------------*/
/*-----------------------  QUATERNION — Coord  -------------------*/


// 将矢量（或三维空间的一点）表示成四元组
QUATERNION VectorToQuaternion(const Coord V);

// 将四元组的虚部用矢量表示
Coord QuaternionToVector(const QUATERNION Q);

// V 经过 Q 变换后为 V', 函数返回 V'
Coord VectorTransform(const Coord V, const QUATERNION Q);
Coord VectorTransform(const Coord V, const ROTATION Q);

#endif