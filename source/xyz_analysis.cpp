/* functions for analysis of xyz_crd class must #include "xyz_analysis.h" */

#include <iostream>
#include <math.h>
using namespace std;

#include "load_crdf.h"
#include "class_log.h" //use of Log() function
#include "xyz_analysis.h"


double xyz_dihedral(int A1,int A2, int A3, int A4, XYZ_CRD * xyz){ 
//generates dihedral angles from 4 cartesian coordinates given atom indexes and xyz crd pointer
   double B1[4],B2[4],B3[4];  //bond vectors                                    B1   B2   B3
   double X_1,Y_1,DIHED,DIHED_DEG; //other vectors  A1---A2---A3---A4  
   //c1[4],c2[4],c3[4],c4[4],
   // [1]==x, [2]==y, [3]==z 
   B1[1]= (xyz->x(A2) - xyz->x(A1) ); 
   B1[2]= (xyz->y(A2) - xyz->y(A1) ); 
   B1[3]= (xyz->z(A2) - xyz->z(A1) ); 
   
   B2[1]= (xyz->x(A3) - xyz->x(A2) ); 
   B2[2]= (xyz->y(A3) - xyz->y(A2) ); 
   B2[3]= (xyz->z(A3) - xyz->z(A2) );   
   
   B3[1]= (xyz->x(A4) - xyz->x(A3) ); 
   B3[2]= (xyz->y(A4) - xyz->y(A3) ); 
   B3[3]= (xyz->z(A4) - xyz->z(A3) ); 
 
   B1[0]=sqrt((B1[1]*B1[1])+(B1[2]*B1[2])+(B1[3]*B1[3]));
   B2[0]=sqrt((B2[1]*B2[1])+(B2[2]*B2[2])+(B2[3]*B2[3]));
   B3[0]=sqrt((B3[1]*B3[1])+(B3[2]*B3[2])+(B3[3]*B3[3]));
 
   //Y_1 = |b2|b1 . [b2 x b3] ; let c1=[b2 x b3]
   double * c1 = cross(B2,B3);
   Y_1= B2[0]*( (B1[1]*c1[1]) + (B1[2]*c1[2]) + (B1[3]*c1[3]) );
   
   //X_1 = [b1 x b2] . [b2 x b3]; let c2 = [b1 x b2]
   double * c2 = cross(B1,B2);
   X_1= ( (c1[1]*c2[1]) + (c1[2]*c2[2]) + (c1[3]*c2[3]) );
   DIHED = atan2(Y_1,X_1);
   //cout << n1u[1] << "  " << n1u[2] << "  "<< n1u[3] << "  " << n2u[1] << "  "<< n2u[2] << "  " << n2u[3] << "  " << m1[0] << endl;
   DIHED_DEG=DIHED*180/3.141592;
   Log(3) << "read atoms "<< A1<<","<<A2<<","<<A3<<","<<A4<<","<< DIHED << " rad, "<< DIHED_DEG<< " deg \n";
   delete c1;
   delete c2;
return DIHED_DEG;
}

double xyzvec_dihedral(int A1,int A2, int A3, int A4, XYZ_VEC_CRD * xyz){
//generates dihedral angles from 4 cartesian coordinates given atom indexes and xyz crd pointer
   double B1[4],B2[4],B3[4];  //bond vectors                                    B1   B2   B3
   double X_1,Y_1,DIHED,DIHED_DEG; //other vectors  A1---A2---A3---A4  
   //c1[4],c2[4],c3[4],c4[4],
   // [1]==x, [2]==y, [3]==z 
   double C[4][3]; //atom coordinates
   for(int an=0; an<4;an++){
	int ai;
	if(an==0) ai = A1; if(an==1) ai = A2;
	if(an==2) ai = A3; if(an==3) ai = A4;
			 if(ai== 0){ C[an][0]=0; C[an][1]=0; C[an][2]=0; } 
		else if(ai==-1){ C[an][0]=1; C[an][1]=0; C[an][2]=0; }
		else if(ai==-2){ C[an][0]=0; C[an][1]=1; C[an][2]=0; }
		else if(ai==-3){ C[an][0]=0; C[an][1]=0; C[an][2]=1; }
		else{
			C[an][0] = xyz->x(ai);
			C[an][1] = xyz->y(ai);
			C[an][2] = xyz->z(ai);
		}
	}
   B1[1]= (C[1][0] - C[0][0] ); 
   B1[2]= (C[1][1] - C[0][1] ); 
   B1[3]= (C[1][2] - C[0][2] ); 
   
   B2[1]= (C[2][0] - C[1][0] ); 
   B2[2]= (C[2][1] - C[1][1] ); 
   B2[3]= (C[2][2] - C[1][2] );   
   
   B3[1]= (C[3][0] - C[2][0] ); 
   B3[2]= (C[3][1] - C[2][1] ); 
   B3[3]= (C[3][2] - C[2][2] ); 
 
   B1[0]=sqrt((B1[1]*B1[1])+(B1[2]*B1[2])+(B1[3]*B1[3]));
   B2[0]=sqrt((B2[1]*B2[1])+(B2[2]*B2[2])+(B2[3]*B2[3]));
   B3[0]=sqrt((B3[1]*B3[1])+(B3[2]*B3[2])+(B3[3]*B3[3]));
 
   //Y_1 = |b2|b1 . [b2 x b3] ; let c1=[b2 x b3]
   double * c1 = cross(B2,B3);
   Y_1= B2[0]*( (B1[1]*c1[1]) + (B1[2]*c1[2]) + (B1[3]*c1[3]) );
   
   //X_1 = [b1 x b2] . [b2 x b3]; let c2 = [b1 x b2]
   double * c2 = cross(B1,B2);
   X_1= ( (c1[1]*c2[1]) + (c1[2]*c2[2]) + (c1[3]*c2[3]) );
   DIHED = atan2(Y_1,X_1);
   //cout << n1u[1] << "  " << n1u[2] << "  "<< n1u[3] << "  " << n2u[1] << "  "<< n2u[2] << "  " << n2u[3] << "  " << m1[0] << endl;
   DIHED_DEG=DIHED*180/3.141592;
   Log(3) << "read atoms "<< A1<<","<<A2<<","<<A3<<","<<A4<<","<< DIHED << " rad, "<< DIHED_DEG<< " deg \n";
   delete c1;
   delete c2;
return DIHED_DEG;
}

double xyzvec_Bangle(int A1,int A2, int A3, XYZ_VEC_CRD * xyz){
//  diagram 1->2->3 p2 at vertex and p21 is len between 2 and 1
//angle calc by arccos((p12^2 +p23^2 -p13^2 )/(2*p21*p23))
   double C[3][3]; //atom coordinates
   for(int an=0; an<3;an++){
	int ai;
	if(an==0) ai = A1; if(an==1) ai = A2;
	if(an==2) ai = A3; 
			 if(ai== 0){ C[an][0]=0; C[an][1]=0; C[an][2]=0; } 
		else if(ai==-1){ C[an][0]=1; C[an][1]=0; C[an][2]=0; }
		else if(ai==-2){ C[an][0]=0; C[an][1]=1; C[an][2]=0; }
		else if(ai==-3){ C[an][0]=0; C[an][1]=0; C[an][2]=1; }
		else{
			C[an][0] = xyz->x(ai);
			C[an][1] = xyz->y(ai);
			C[an][2] = xyz->z(ai);
		}
	}
	double B1[4],B2[4],B3[4]; 
    B1[1]= (C[1][0] - C[0][0] ); 
    B1[2]= (C[1][1] - C[0][1] ); 
    B1[3]= (C[1][2] - C[0][2] ); 
    
    B2[1]= (C[2][0] - C[1][0] ); 
    B2[2]= (C[2][1] - C[1][1] ); 
    B2[3]= (C[2][2] - C[1][2] );   
    
    B3[1]= (C[0][0] - C[2][0] ); 
    B3[2]= (C[0][1] - C[2][1] ); 
    B3[3]= (C[0][2] - C[2][2] ); 
    
    B1[0]=sqrt((B1[1]*B1[1])+(B1[2]*B1[2])+(B1[3]*B1[3]));
    B2[0]=sqrt((B2[1]*B2[1])+(B2[2]*B2[2])+(B2[3]*B2[3]));
    B3[0]=sqrt((B3[1]*B3[1])+(B3[2]*B3[2])+(B3[3]*B3[3]));
	double ANG = (B1[0]*B1[0]) + (B2[0]*B2[0]) - (B3[0]*B3[0]);
	ANG = ANG / (2 * B1[0] * B2[0]) ;
	ANG = acos(ANG);
return ANG*180.0/3.14159265359;
}

double xyzvec_dist(int A1, int A2,  XYZ_VEC_CRD * xyz){
	//return distance between 2 atoms
	double dx = xyz->x(A1) - xyz->x(A2); 
	double dy = xyz->y(A1) - xyz->y(A2);
	double dz = xyz->z(A1) - xyz->z(A2);
	double D=sqrt((dx*dx)+(dy*dy)+(dz*dz));
return D;
}


 double * cross(double V_1[], double V_2[]){
    double * V_3 = new double[4];
    V_3[1]= (V_1[2]*V_2[3]) - (V_2[2]*V_1[3]);
	V_3[2]=-(V_1[1]*V_2[3]) + (V_2[1]*V_1[3]);
	V_3[3]= (V_1[1]*V_2[2]) - (V_2[1]*V_1[2]);
    V_3[0] = sqrt( (V_3[1]*V_3[1])+(V_3[2]*V_3[2]) + (V_3[3]*V_3[3]) );
return V_3;
 }
 
 double xyz_checkRange(XYZ_VEC_CRD * crd, double X, double Y, double Z){
	//return the minimum distance between xyz point and ref structure
	double Dmin=9e99;
	for(int jj=1; jj<= crd->n(); jj++){
		double dx=crd->x(jj)-X;
		double dy=crd->y(jj)-Y;
		double dz=crd->z(jj)-Z;
		double D=sqrt((dx*dx)+(dy*dy)+(dz*dz));
		if(D < Dmin) Dmin = D;
	}
	Log(3)<<"F:xyz_checkRange: Min Dist "<<Dmin<<" A\n";
 return Dmin;
 }