#ifndef HJ_TBLAS_H_
#define HJ_TBLAS_H_

namespace hj { namespace tblas {

template <typename V0, typename V1, typename V2>
inline void cross(const V0 &v0, const V1 &v1, V2 &v2) {
	v2[0] = v0[1]*v1[2]-v0[2]*v1[1];
	v2[1] = v0[2]*v1[0]-v0[0]*v1[2];
	v2[2] = v0[0]*v1[1]-v0[1]*v1[0];
}

template <typename T, int size, typename V1, typename V2>
inline T dot(const V1 &v1, const V2 &v2) {
	T rtn(0);
	for(int i = 0; i < size; ++i) rtn += v1[i]*v2[i];
	return rtn;
}

template <typename T, int size, typename V>
inline T norm(const V &v) {
	return sqrt(dot<T, size>(v, v));
}

template <typename T>
inline int sign(T v) {
	if(v > 0) return 1;
	if(v < 0) return -1;
	return 0;
}

}}

#endif
