#ifndef XYZ_ANALYSIS_H_
#define XYZ_ANALYSIS_H_

double xyz_dihedral(int A1,int A2, int A3, int A4, XYZ_CRD * xyz); //generates dihedral angles from 4 cartesian coordinates given atom indexes and xyz crd pointer
double xyzvec_dihedral(int A1,int A2, int A3, int A4, XYZ_VEC_CRD * xyz);
double xyzvec_Bangle(int A1,int A2, int A3, XYZ_VEC_CRD * xyz);
double xyzvec_dist(int A1, int A2,  XYZ_VEC_CRD * xyz);

double * cross(double V_1[], double V_2[]); // do cross product of doble array with 4 components [0]=magnitude [1]=x-comp [2]=y-comp [3]=z-comp

 double xyz_checkRange(XYZ_VEC_CRD * crd, double X, double Y, double Z);

#endif