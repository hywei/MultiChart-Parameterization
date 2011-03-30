#ifndef HJ_QUATERNION_H_
#define HJ_QUATERNION_H_

#include <hjlib/math/tblas.h>
#include <cmath>
#include <limits>
#include <algorithm>

namespace hj { namespace quaternion {

// Q is [w x y z]
template <typename Q>
inline const Q &eye(Q &q) {
	q[0] = 1; q[1] = 0; q[2] = 0; q[3] = 0;
	return q;
}

// output is stored in q
template <typename T, typename Q>
inline bool axis_angle(Q &q) {
	const T n = hj::tblas::norm<T, 4>(q);
	if(n < 1e-5) return false;
	for(int i = 0; i < 4; ++i) q[i] /= n;

    const T cos_a = q[0];
    q[0] = acos(cos_a)*2;
    T sin_a = sqrt( 1.0 - cos_a*cos_a );
    if(fabs(sin_a) < 0.0005)
		sin_a = 1;
    q[1] /= sin_a;
    q[2] /= sin_a;
    q[3] /= sin_a;
	return true;
}

template <typename Q, typename A, typename T>
inline const Q& axis_angle(const A &axis, T theta, Q &q) {
	T sin_half_theta = sin(theta/2);
	assert(sin_half_theta*sin_half_theta <= 1.0f);
	q[0] = sqrt(1-sin_half_theta*sin_half_theta);
	q[1] = axis[0]*sin_half_theta,
	q[2] = axis[1]*sin_half_theta,
	q[3] = axis[2]*sin_half_theta;
	return q;
}

// minimal q rotate v0 to v1
// NOTE: v0, v1 must be normalized
template <typename T, typename Q, typename V1, typename V2>
inline const Q& min_rot(const V1 &v0, const V2 &v1, Q &q)
{
	const T eps = std::numeric_limits<T>::epsilon();
	T axis[3];
	hj::tblas::cross(v0, v1, axis);
	const T len_axis = hj::tblas::norm<T, 3>(axis);
	T cos_theta = hj::tblas::dot<T, 3>(v0, v1), sin_theta = len_axis;
	if(len_axis < eps) {	// parallel
		q[1] = q[2] = q[3] = 0;
		if(cos_theta > 0)
			q[0] = 1;
		else
			q[0] = 0;
		return q;
	}
	axis[0] /= len_axis;
	axis[1] /= len_axis;
	axis[2] /= len_axis;

	if(cos_theta > 1) cos_theta = 1;
	if(cos_theta <-1) cos_theta = -1;
	const T cos_half_theta = sqrt((1+cos_theta)/2);
	const T sin_half_theta = sqrt(1-cos_half_theta*cos_half_theta)*hj::tblas::sign(sin_theta);

	q[0] = cos_half_theta;
	q[1] = axis[0]*sin_half_theta,
	q[2] = axis[1]*sin_half_theta,
	q[3] = axis[2]*sin_half_theta;
	return q;
}

// Q is [w x y z]
// M is pre-alloc column major 3x3 matrix layout object
template <typename T, typename Q, typename M>
const M& quat2m33(const Q &q, M &m)
{
	typedef T value_type;

    const value_type a = q[0], b = q[1], c = q[2], d = q[3];    
    const value_type aa = a*a, ab = a*b, ac = a*c, ad = a*d,
		bb = b*b, bc = b*c, bd = b*d,
		cc = c*c, cd = c*d,
		dd = d*d;
    const value_type norme_carre = aa+bb+cc+dd;
    
	if(norme_carre <= std::numeric_limits<value_type>::epsilon()) {
		m[1] = m[2] = m[3] = m[5] = m[6] = m[7] = 0;
		m[0] = m[4] = m[8] = -1;
	}
	else {
		m[0] = (aa+bb-cc-dd)/norme_carre;
		m[1] = 2*(ad+bc)/norme_carre;
		m[2] = 2*(-ac+bd)/norme_carre;
		m[3] = 2*(-ad+bc)/norme_carre;
		m[4] = (aa-bb+cc-dd)/norme_carre;
		m[5] = 2*(ab+cd)/norme_carre;
		m[6] = 2*(ac+bd)/norme_carre;
		m[7] = 2*(-ab+cd)/norme_carre;
		m[8] = (aa-bb-cc+dd)/norme_carre;
	}
	return m;
}

// Q is [w x y z]
// M is pre-alloc column major 3x3 rotation matrix layout array
template<typename T, typename Q, typename M>
const Q &m332quat(const M& m, Q &q)
{
	typedef T value_type;
	const value_type trace[4] = {
		+m[0]+m[4]+m[8],
		+m[0]-m[4]-m[8],
		-m[0]+m[4]-m[8],
		-m[0]-m[4]+m[8]
	};
	const int max_idx = static_cast<int>(std::max_element(trace, trace+4)-trace);
	const value_type s = sqrt(trace[max_idx]+1)*2;
	switch(max_idx) {
		case 0:
			q[0] = s/4;
			q[1] = (m[5]-m[7])/s;
			q[2] = (m[6]-m[2])/s;
			q[3] = (m[1]-m[3])/s;
			break;
		case 1:
			q[0] = (m[5]-m[7])/s;
			q[1] = s/4;
			q[2] = (m[3]+m[1])/s;
			q[3] = (m[2]+m[6])/s;
			break;
		case 2:
			q[0] = (m[6]-m[2])/s;
			q[1] = (m[3]+m[1])/s;
			q[2] = s/4;
			q[3] = (m[7]+m[5])/s;
			break;
		case 3:
			q[0] = (m[1]-m[3])/s;
			q[1] = (m[2]+m[6])/s;
			q[2] = (m[7]+m[5])/s;
			q[3] = s/4;
			break;
	}
    if (q[0] <0) {
        q[0] = -q[0];
        q[1] = -q[1];
        q[2] = -q[2];
        q[3] = -q[3];
    }
	return q;
}

}}

#endif
