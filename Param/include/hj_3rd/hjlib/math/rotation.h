#ifndef HJ_ROTATION_H_
#define HJ_ROTATION_H_

#include "conf.h"

#ifdef __cplusplus
extern "C" {
#endif

// axis is a normalized length 3 vector
// rot is column major 3x3 matrix
// output should be pre-allocated
void HJ_MATH_API axis_angle_to_rot(double *axis, double angle, double *rot);

/**
   @param axisangle will be normalized

   @return 1 means too small rotation (angle*angle < eps)
 */
int HJ_MATH_API axisangle_to_rot(double *axisangle, double *rot, double eps);

/**
   @param rot will not be modified
 */
void HJ_MATH_API rot_to_axis_angle_Rodrigiues(double *rot, double *axis, double *angle, double eps);

void HJ_MATH_API rot_to_axis_angle_direct(double *rot, double *axis, double *angle, double eps);

#ifdef __cplusplus
}
#endif

#endif
