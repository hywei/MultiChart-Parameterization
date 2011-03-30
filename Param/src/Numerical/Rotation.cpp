/*_____________________________________________________
 |
 |	Intra3D Layer1 :  Algebra 构件组成之一
 |                                          
 |  文件：Rotation.cpp
 |
 |  功能：四元组运算
 |
 |  开发：林锐 ，1999/01/10
 |
 |	源程序测试：进行了单步跟踪测试
 |
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


#include <math.h>
#include <assert.h>
#include "Rotation.h"


/*_____________________________________________________
 |                                                     
 |          Quaternion & Rotation Functions     
 |
 |        已进行逐步跟踪测试，林锐，1999/01/16
 |_____________________________________________________*/

/*_____________________________________________________________
 |                                                     
 | 有关 Quaternion 与 Rotation 的原理请参考：     
 |
 | (1) Laura Downs, Using Quaternions to Represent Rotation,
 |     CS184 at UC Berkeley
 |     http://www.cs.berkeley.edu/~laura			
 |
 | (2) Nick Bobick, Rotating Objects Using Quaternions,
 |     Game Developer, July 3, 1998
 |     http://www.gamasutra.com/features/programming/19980703			
 |_____________________________________________________________*/

// 旋转误差
//const double DELTA_ROT=0.000001f;
const double DELTA_ROT=1.0E-10;

// ROTATION 相乘，返回值 = R1 * R2
// 旋转变换顺序：先执行 R1 旋转，后执行 R2 旋转
// 要求 R1 R2 先归一化
ROTATION RotationMultiply( ROTATION R1,  ROTATION R2)
{
	// 原始角度为 0
	if( (R1.angle>-DELTA_ROT) && (R1.angle<DELTA_ROT) )
	{
		return R2;
	}
	else if( (R2.angle>-DELTA_ROT) && (R2.angle<DELTA_ROT) )
	{
		return R1;
	}
	
	// 原始角度不为 0
	ROTATION   R; 
	QUATERNION Q1=RotationToQuaternion(R1);
	QUATERNION Q2=RotationToQuaternion(R2);
	QUATERNION Q = Q2*Q1; // 注意顺序
	R=QuaternionToRotation(Q);
	return R;
}

ROTATION operator * ( ROTATION R1,  ROTATION R2)
{
	return RotationMultiply(R1, R2);
}


// V 经过 R 变换后为 V', 函数返回 V'
// 要求 R 先归一化
Coord VectorTransform(const Coord V, const ROTATION R)
{
	if( (R.angle>-DELTA_ROT) && (R.angle<DELTA_ROT) ) return V;

	QUATERNION Q;
	Q=RotationToQuaternion(R);
	return VectorTransform(V, Q);
}

// 将 ROTATION 结构表示成 Matrix
// 要求 R 先归一化
void RotationToMatrix(double M[16], const ROTATION R)
{
	if( (R.angle>-DELTA_ROT) && (R.angle<DELTA_ROT) ) 
	{
		M[0] = M[5] = M[10] = M[15] = 1.0;
		M[1] = M[2] = M[3] = M[4] = M[6] = M[7] = M[8] = M[9] = M[11] = M[12] = M[13] = M[14] = 0.0;
		
		return ;
	}
	else
		QuaternionToMatrix(M, RotationToQuaternion(R));
}

/*_______________ 不想了解四元组细节的程序员不必往下看 _____________________*/




// 四元组求模，返回值 = |A|
double QuaternionMagnitude(const QUATERNION A)
{
	return double(sqrt(A.x*A.x + A.y*A.y + A.z*A.z +A.w*A.w));
}

// 四元组归一化
// 如果 |A|=0，输出值 = A；输出值 = A/(|A|)
void  QuaternionNormalize(QUATERNION *A)
{
	double magnitude = double(sqrt(A->x*A->x + A->y*A->y + A->z*A->z +A->w*A->w));

	if(magnitude >= DELTA_ROT)
    {
	  A->x = A->x/magnitude;
	  A->y = A->y/magnitude;
	  A->z = A->z/magnitude;
	  A->w = A->w/magnitude;
    }
}

// 四元组求逆
// 如果 |A|=0，输出值 = A；否则输出值 = A 的逆
void  QuaternionInverse(QUATERNION *A)
{
	double magnitude2 = A->x*A->x + A->y*A->y + A->z*A->z + A->w*A->w;

	if(magnitude2 >= DELTA_ROT)
    {
	  A->x = -A->x/magnitude2;
	  A->y = -A->y/magnitude2;
	  A->z = -A->z/magnitude2;
	  A->w =  A->w/magnitude2;
    }
}

// 四元组共扼，输出值 = A 的共扼
void  QuaternionConjugate(QUATERNION *A)
{
	A->x = -A->x;
	A->y = -A->y;
	A->z = -A->z;
}

// 四元组相加，返回值 = A + B
QUATERNION QuaternionAdd(const QUATERNION A, const QUATERNION B)
{
	QUATERNION Q;
	Q.x = A.x + B.x;
	Q.y = A.y + B.y;
	Q.z = A.z + B.z;
	Q.w = A.w + B.w;
	return Q;
}

QUATERNION operator+(const QUATERNION A, const QUATERNION B)
{
	QUATERNION Q;
	Q.x = A.x + B.x;
	Q.y = A.y + B.y;
	Q.z = A.z + B.z;
	Q.w = A.w + B.w;
	return Q;
}

// 四元组相减，返回值 = A - B
QUATERNION QuaternionSub(const QUATERNION A, const QUATERNION B)
{
	QUATERNION Q;
	Q.x = A.x - B.x;
	Q.y = A.y - B.y;
	Q.z = A.z - B.z;
	Q.w = A.w - B.w;
	return Q;
}

QUATERNION operator-(const QUATERNION A, const QUATERNION B)
{
	QUATERNION Q;
	Q.x = A.x - B.x;
	Q.y = A.y - B.y;
	Q.z = A.z - B.z;
	Q.w = A.w - B.w;
	return Q;
}


// 四元组缩放，返回值 = s * A 
QUATERNION QuaternionScale(const QUATERNION A, const double s)
{
	QUATERNION Q;
	Q.x = s * A.x;
	Q.y = s * A.y;
	Q.z = s * A.z;
	Q.w = s * A.w;
	return Q;
}

QUATERNION operator * (const QUATERNION A, const double s)
{
	QUATERNION Q;
	Q.x = s * A.x;
	Q.y = s * A.y;
	Q.z = s * A.z;
	Q.w = s * A.w;
	return Q;
}

QUATERNION QuaternionScale(const double s, const QUATERNION A)
{
	QUATERNION Q;
	Q.x = s * A.x;
	Q.y = s * A.y;
	Q.z = s * A.z;
	Q.w = s * A.w;
	return Q;
}

QUATERNION operator * (const double s, const QUATERNION A)
{
	QUATERNION Q;
	Q.x = s * A.x;
	Q.y = s * A.y;
	Q.z = s * A.z;
	Q.w = s * A.w;
	return Q;
}

// 四元组相乘，返回值 = q1 * q2
QUATERNION QuaternionMultiply(const QUATERNION q1, const QUATERNION q2)
{
// 林锐编写，已测试正确
	QUATERNION Q;
    Q.x =	q1.w * q2.x + q1.x * q2.w +	q1.y * q2.z - q1.z * q2.y;
    Q.y =	q1.w * q2.y + q1.y * q2.w +	q1.z * q2.x - q1.x * q2.z;
    Q.z =	q1.w * q2.z + q1.z * q2.w +	q1.x * q2.y - q1.y * q2.x;
    Q.w =	q1.w * q2.w - q1.x * q2.x -	q1.y * q2.y - q1.z * q2.z;
	return Q;
}

// 四元组相乘，返回值 = q1 * q2
QUATERNION operator * (const QUATERNION q1, const QUATERNION q2)
{
// 林锐编写，已测试正确
	QUATERNION Q;
    Q.x =	q1.w * q2.x + q1.x * q2.w +	q1.y * q2.z - q1.z * q2.y;
    Q.y =	q1.w * q2.y + q1.y * q2.w +	q1.z * q2.x - q1.x * q2.z;
    Q.z =	q1.w * q2.z + q1.z * q2.w +	q1.x * q2.y - q1.y * q2.x;
    Q.w =	q1.w * q2.w - q1.x * q2.x -	q1.y * q2.y - q1.z * q2.z;
	return Q;
/*
// 林锐编写，已测试正确
	double  s, s1, s2;
	Coord v, v1, v2;

	s1 = q1.w;
	s2 = q2.w;
	v1.x=q1.x; v1.y=q1.y; v1.z=q1.z;
	v2.x=q2.x; v2.y=q2.y; v2.z=q2.z;

	s  = s1*s2 - v1^v2;
	v  = s1*v2 + s2*v1 + v1*v2;
	QUATERNION Q;
	Q.w= s;
	Q.x= v.x;
	Q.y= v.y;
	Q.z= v.z;
	return Q;
*/
/*
// Jeff Lander 编写，变成了 Q2*Q1
	QUATERNION Q;

    Q.x =	q2.w * q1.x + q2.x * q1.w +
			q2.y * q1.z - q2.z * q1.y;
    Q.y =	q2.w * q1.y + q2.y * q1.w +
			q2.z * q1.x - q2.x * q1.z;
    Q.z =	q2.w * q1.z + q2.z * q1.w +
			q2.x * q1.y - q2.y * q1.x;
    Q.w =	q2.w * q1.w - q2.x * q1.x -
			q2.y * q1.y - q2.z * q1.z;
	return Q;
*/

/*
// Nick Bobick 编写，变成了 Q2*Q1

	double A, B, C, D, E, F, G, H;

	A = (q1.w + q1.x) * (q2.w + q2.x);
	B = (q1.z - q1.y) * (q2.y - q2.z);
	C = (q1.x - q1.w) * (q2.y - q2.z);
	D = (q1.y + q1.z) * (q2.x - q2.w);
	E = (q1.x + q1.z) * (q2.x + q2.y);
	F = (q1.x - q1.z) * (q2.x - q2.y);
	G = (q1.w + q1.y) * (q2.w - q2.z);
	H = (q1.w - q1.y) * (q2.w + q2.z);

	QUATERNION Q;
	Q.w =  B + (-E - F + G + H)/2.0f;
	Q.x =  A - ( E + F + G + H)/2.0f; 
	Q.y = -C + ( E - F + G - H)/2.0f;
	Q.z = -D + ( E - F - G + H)/2.0f;
	return Q;
*/
}

// Nick Bobick 编写
// Purpose:		Spherical Linear Interpolation Between two Quaternions
// Arguments:	Two Quaternions, blend factor
QUATERNION QuaternionSlerp(const QUATERNION from, const QUATERNION to, double t)
{
		double  to1[4];
        double  omega, cosom, sinom, scale0, scale1;

        // calc cosine
        cosom = from.x * to.x + from.y * to.y + from.z * to.z
			       + from.w * to.w;

        // adjust signs (if necessary)
        if ( cosom <0.0 ){ cosom = -cosom; to1[0] = - to.x;
		to1[1] = - to.y;
		to1[2] = - to.z;
		to1[3] = - to.w;
        } else  {
		to1[0] = to.x;
		to1[1] = to.y;
		to1[2] = to.z;
		to1[3] = to.w;
        }

        // calculate coefficients

       if ( (1.0 - cosom) > DELTA_ROT ) {
                // standard case (slerp)
                omega = double(acos(cosom));
                sinom = double(sin(omega));
                scale0 = double(sin((1.0 - t) * omega) / sinom);
                scale1 = double(sin(t * omega) / sinom);

        } else {        
    // "from" and "to" quaternions are very close 
	    //  ... so we can do a linear interpolation
                scale0 = 1.0f - t;
                scale1 = t;
        }
	// calculate final values
	QUATERNION res;
	res.x = scale0 * from.x + scale1 * to1[0];
	res.y = scale0 * from.y + scale1 * to1[1];
	res.z = scale0 * from.z + scale1 * to1[2];
	res.w = scale0 * from.w + scale1 * to1[3];
	return res;
}


/*_____________________________________________________
 |                                                     
 |   为提高计算性能，以下变换函数均假定
 |
 |   输入参数 R Q V 已经进行了 Normalize 处理    
 |_____________________________________________________*/



/*-----------------------------------------------*/
/*-----------  QUATERNION — ROTATION  ----------*/

// 将 ROTATION 结构表示成 QUATERNION
// 假定 ROTATION 参数已经进行了归一化处理
QUATERNION RotationToQuaternion(const ROTATION R)
{
	QUATERNION Q;
	double cosValue, sinValue, theta;

	theta    = R.angle/180.0f*3.14159f;
	cosValue = double(cos(theta/2.0f));  
	sinValue = double(sin(theta/2.0f));  

	Q.x = sinValue*R.axis[0];
	Q.y = sinValue*R.axis[1];
	Q.z = sinValue*R.axis[2];
	Q.w = cosValue;

	return Q;
}

// 将 QUATERNION 结构表示成 ROTATION
// 假定 QUATERNION 参数已经进行了归一化处理
ROTATION   QuaternionToRotation(const QUATERNION Q)
{
	ROTATION   R; 
	if( ((Q.w>1-DELTA_ROT)&&(Q.w<1+DELTA_ROT)) || 
		((Q.w>-1-DELTA_ROT)&&(Q.w<-1+DELTA_ROT)) ) 
	{
		R.angle=0.0f;
		R.axis[0]=0.0f;
		R.axis[1]=0.0f;
		R.axis[2]=1.0f;
		return R;
	}

	double sinValue, halfTheta;
	halfTheta= double(acos(Q.w));
	sinValue = double(sin(halfTheta));  
	
//	assert(sinValue >= DELTA_ROT);

	if(sinValue <= DELTA_ROT)
	{
		R.angle=0.0f;
		R.axis[0]=0.0f;
		R.axis[1]=0.0f;
		R.axis[2]=1.0f;
		return R;
	}

	R.angle  = halfTheta * 2.0f * 180.0f /3.14159f;
	R.axis[0] = Q.x/sinValue;
	R.axis[1] = Q.y/sinValue;
	R.axis[2] = Q.z/sinValue;
//	VectorNormalize(R.axis);  // R.axis 必是单位矢量
	return R;
}


/*-------------------------------------------------------------*/
/*-------------------  QUATERNION — Matrix  ------------------*/

// Nick Bobick 编写
// 将 Matrix 结构表示成 QUATERNION
// 要求 Matrix 是一种旋转矩阵,否则不能得到正确结果
void  MatrixToQuaternion(QUATERNION *quat, const double M[16])
{
	double tr, s, q[4];
	int   i, j, k;
	int	  nxt[3] = {1, 2, 0};

	double m[4][4];
	for(i=0; i<4; i++)
	for(j=0; j<4; j++)
		m[i][j]=M[i*4+j];	

	tr = m[0][0] + m[1][1] + m[2][2];

	// check the diagonal
	if (tr > 0.0) 
	{
		s = double(sqrt (tr + 1.0));   // s 为什么 !=  -sqrt (tr + 1.0)
		quat->w = s / 2.0f;
		s = 0.5f / s;
		quat->x = (m[1][2] - m[2][1]) * s;
		quat->y = (m[2][0] - m[0][2]) * s;
		quat->z = (m[0][1] - m[1][0]) * s;
	} 
	else 
	{		
	 // diagonal is negative
    	i = 0;
        if (m[1][1] > m[0][0]) i = 1;
	    if (m[2][2] > m[i][i]) i = 2;
        j = nxt[i];
        k = nxt[j];
        s = double(sqrt ((m[i][i] - (m[j][j] + m[k][k])) + 1.0f));
	    q[i] = s * 0.5f;
        if (s != 0.0f) s = 0.5f / s;
	    q[3] = (m[j][k] - m[k][j]) * s;
        q[j] = (m[i][j] + m[j][i]) * s;
        q[k] = (m[i][k] + m[k][i]) * s;
		quat->x = q[0];
		quat->y = q[1];
		quat->z = q[2];
		quat->w = q[3];
	}
}

void MatrixToQuaternion2(QUATERNION& quat, const double R[3][3])
{
	double trace = 1 + R[0][0] + R[1][1] + R[2][2];
	double s;
	if(trace > 1.0e-8)
	{
		s = sqrt(trace)*2.0;
		quat.x = ( R[2][1] - R[1][2] ) / s;
		quat.y = ( R[0][2] - R[2][0] ) / s;
		quat.z = ( R[1][0] - R[0][1] ) / s;
		quat.w = 0.25 * s;
	}
	else
	{
		if ( R[0][0] > R[1][1] && R[0][0] > R[2][2] )  {	// Column 0: 
			s  = sqrt( 1.0 + R[0][0] - R[1][1] - R[2][2] ) * 2;
			quat.x = 0.25 * s;
			quat.y = (R[1][0] + R[0][1] ) / s;
			quat.z = (R[0][2] + R[2][0] ) / s;
			quat.w = (R[2][1] - R[1][2] ) / s;
		} else if ( R[1][1] > R[2][2] ) {			// Column 1: 
			s  = sqrt( 1.0 + R[1][1] - R[0][0] - R[2][2] ) * 2;
			quat.x = (R[1][0] + R[0][1] ) / s;
			quat.y = 0.25 * s;
			quat.z = (R[2][1] + R[1][2] ) / s;
			quat.w = (R[0][2] - R[2][0] ) / s;
		} else {						// Column 2:
			s  = sqrt( 1.0 + R[2][2] - R[0][0] - R[1][1] ) * 2;
			quat.x = (R[0][2] + R[2][0] ) / s;
			quat.y = (R[2][1] + R[1][2] ) / s;
			quat.z = 0.25 * s;
			quat.w = (R[1][0] - R[0][1] ) / s;
		}
	}
//	if(trace > 1.0e-8)
//	{
//		s = sqrt(trace)*2.0;
//		quat.x = ( R[1][2] - R[2][1] ) / s;
//		quat.y = ( R[2][0] - R[0][2] ) / s;
//		quat.z = ( R[0][1] - R[1][0] ) / s;
//		quat.w = 0.25 * s;
//	}
//	else
//	{
//		if ( R[0][0] > R[1][1] && R[0][0] > R[2][2] )  {	// Column 0: 
//			s  = sqrt( 1.0 + R[0][0] - R[1][1] - R[2][2] ) * 2;
//			quat.x = 0.25 * s;
//			quat.y = (R[0][1] + R[1][0] ) / s;
//			quat.z = (R[2][0] + R[0][2] ) / s;
//			quat.w = (R[1][2] - R[2][1] ) / s;
//		} else if ( R[1][1] > R[2][2] ) {			// Column 1: 
//			s  = sqrt( 1.0 + R[1][1] - R[0][0] - R[2][2] ) * 2;
//			quat.x = (R[0][1] + R[1][0] ) / s;
//			quat.y = 0.25 * s;
//			quat.z = (R[1][2] + R[2][1] ) / s;
//			quat.w = (R[2][0] - R[0][2] ) / s;
//		} else {						// Column 2:
//			s  = sqrt( 1.0 + R[2][2] - R[0][0] - R[1][1] ) * 2;
//			quat.x = (R[2][0] + R[0][2] ) / s;
//			quat.y = (R[1][2] + R[2][1] ) / s;
//			quat.z = 0.25 * s;
//			quat.w = (R[0][1] - R[1][0] ) / s;
//		}
//	}
}

void MatrixToQuaternion(QUATERNION *quat, const double R[3][3])
{
	double M[16];
	int i,j;
	for(i = 0; i < 3; i ++)
		for(j = 0; j < 3; j ++)
			M[i*4+j] = R[i][j];
	M[3] = M[1*4+3] = M[2*4+3] = M[3*4+0] = M[3*4+1] = M[3*4+2] = 0.0;
	M[3*4+3] = 1.0;
	MatrixToQuaternion(quat,M);
}

// 将 QUATERNION 结构表示成 Matrix
void  QuaternionToMatrix(double M[16], const QUATERNION quat)
{
  double m[4][4];  

  double wx, wy, wz, xx, yy, yz, xy, xz, zz, x2, y2, z2;

  // calculate coefficients
  x2 = quat.x * 2.0f; 
  y2 = quat.y * 2.0f; 
  z2 = quat.z * 2.0f;
  xx = quat.x * x2;   xy = quat.x * y2;   xz = quat.x * z2;
  yy = quat.y * y2;   yz = quat.y * z2;   zz = quat.z * z2;
  wx = quat.w * x2;   wy = quat.w * y2;   wz = quat.w * z2;

  m[0][0] = 1.0f - (yy + zz); 	
  m[0][1] = xy + wz;
  m[0][2] = xz - wy;		
  m[0][3] = 0.0f;
 
  m[1][0] = xy - wz;		
  m[1][1] = 1.0f - (xx + zz);
  m[1][2] = yz + wx;		
  m[1][3] = 0.0f;

  m[2][0] = xz + wy;		
  m[2][1] = yz - wx;
  m[2][2] = 1.0f - (xx + yy);		
  m[2][3] = 0.0f;

  m[3][0] = 0;			
  m[3][1] = 0;
  m[3][2] = 0;			
  m[3][3] = 1;

  for(int i=0; i<4; i++)
  for(int j=0; j<4; j++)
	  M[i*4+j]=m[i][j];	
}

void QuaternionToMatrix(double R[3][3], const QUATERNION quat)
{
	double wx, wy, wz, xx, yy, yz, xy, xz, zz, x2, y2, z2;
	
	// calculate coefficients
	x2 = quat.x * 2.0f; 
	y2 = quat.y * 2.0f; 
	z2 = quat.z * 2.0f;
	xx = quat.x * x2;   xy = quat.x * y2;   xz = quat.x * z2;
	yy = quat.y * y2;   yz = quat.y * z2;   zz = quat.z * z2;
	wx = quat.w * x2;   wy = quat.w * y2;   wz = quat.w * z2;

	// Transpose
//	R[0][0] = 1.0f - (yy + zz);
//	R[1][0] = xy - wz;
//	R[2][0] = xz + wy;
//	
//	R[0][1] = xy + wz;
//	R[1][1] = 1.0f - (xx + zz);
//	R[2][1] = yz - wx;
//	
//	R[0][2] = xz - wy;
//	R[1][2] = yz + wx;
//	R[2][2] = 1.0f - (xx + yy);
	
	R[0][0] = 1.0f - (yy + zz);
	R[0][1] = xy - wz;
	R[0][2] = xz + wy;
	
	R[1][0] = xy + wz;
	R[1][1] = 1.0f - (xx + zz);
	R[1][2] = yz - wx;
	
	R[2][0] = xz - wy;
	R[2][1] = yz + wx;
	R[2][2] = 1.0f - (xx + yy);
}

// 将 QUATERNION 结构表示成 Matrix
// The matrix is represented in a column major format (OpenGL format) 
void  QuaternionToMatrix_OpenGL(double M[16], const QUATERNION quat)
{
  double m[4][4];  //column major format (OpenGL format) 

  double wx, wy, wz, xx, yy, yz, xy, xz, zz, x2, y2, z2;

  // calculate coefficients
  x2 = quat.x * 2.0f; 
  y2 = quat.y * 2.0f; 
  z2 = quat.z * 2.0f;
  xx = quat.x * x2;   xy = quat.x * y2;   xz = quat.x * z2;
  yy = quat.y * y2;   yz = quat.y * z2;   zz = quat.z * z2;
  wx = quat.w * x2;   wy = quat.w * y2;   wz = quat.w * z2;

  m[0][0] = 1.0f - (yy + zz); 	
  m[0][1] = xy - wz;
  m[0][2] = xz + wy;		
  m[0][3] = 0.0f;
 
  m[1][0] = xy + wz;		
  m[1][1] = 1.0f - (xx + zz);
  m[1][2] = yz - wx;		
  m[1][3] = 0.0f;

  m[2][0] = xz - wy;		
  m[2][1] = yz + wx;
  m[2][2] = 1.0f - (xx + yy);		
  m[2][3] = 0.0f;

  m[3][0] = 0;			
  m[3][1] = 0;
  m[3][2] = 0;			
  m[3][3] = 1;

  for(int i=0; i<4; i++)
  for(int j=0; j<4; j++)
	  M[i*4+j]=m[i][j];	//column major format (OpenGL format) 
}

/*-----------------------------------------------*/
/*-----------  QUATERNION — Coord  ------------*/
 
// 将矢量（或三维空间的一点）表示成四元组
QUATERNION VectorToQuaternion(const Coord V)
{
	QUATERNION Q;
    Q.x = V[0];
    Q.y = V[1];
    Q.z = V[2];
    Q.w = 0.0f;
	return Q;
}

// 将四元组的虚部用矢量表示
Coord QuaternionToVector(const QUATERNION Q)
{
	Coord V;
    V[0] = Q.x;
    V[1] = Q.y;
    V[2] = Q.z;
	return V;
}

// V 经过 Q 变换后为 V', 函数返回 V'
Coord VectorTransform(const Coord V, const QUATERNION Q)
{
	if( ((Q.w>1-DELTA_ROT)&&(Q.w<1+DELTA_ROT)) || 
		((Q.w>-1-DELTA_ROT)&&(Q.w<-1+DELTA_ROT)) ) return V;

	QUATERNION A, B, C;
	A=VectorToQuaternion(V);
	B=Q;
	QuaternionInverse(&B);
	C=Q*A*B;
	return QuaternionToVector(C);
}


