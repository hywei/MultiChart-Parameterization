#ifndef TRIANGULARMESH_H_
#define TRIANGULARMESH_H_

#ifdef __cplusplus
extern "C" {
#endif

	struct tri_mesh_3d
	{
		double *vertex;
		int *face;
		int vert_num, face_num;
	};


#ifdef __cplusplus
}
#endif

#endif