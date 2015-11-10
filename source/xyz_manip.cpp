#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <math.h>
#include <vector>
#include <iomanip>
using namespace std;

#include "load_crdf.h"
#include "xyz_manip.h"
#include "class_log.h" //use of Log() function

/*
Functions for manipulating xyz coordinates through the XYZ_CRD class
*/


XYZ_CRD * strip_mol(XYZ_CRD * xyz, int str_ind, int e_ind, int b_ind, double b_len, string a_type){
	// strips away all atoms between str_ind and e_ind (the last atom by defaul) 
	//and replaces it with a 'capping' atom of type 'a_type' with a bond distance of b_len with atom b_ind. 
	XYZ_CRD * ext_crd = new XYZ_CRD;
	int num_init = xyz->n(); 
	int wrtnum=0;
	if(e_ind < 1){ 
		e_ind = num_init;		
	}
	if( e_ind <= str_ind ){
		Log(-1,"strip_tail") <<" improper start and end indecies for extraction. End must be higher than start!\n";
		return ext_crd;
	
	}
	// copy crd up to str_ind
	for( int jj = 1; jj < str_ind; jj++){
		wrtnum++;
		ext_crd->sett(wrtnum, xyz->type(jj) ) ;
		ext_crd->setx(wrtnum, xyz->x(jj) );
		ext_crd->sety(wrtnum, xyz->y(jj) );
		ext_crd->setz(wrtnum, xyz->z(jj) );
		
	}
		// must add in new capping atom bond vector is x_new - x_bond i_hat ...
	//get bond vector for old atom position
	double vec[4]={0,0,0,0};
	vec[1] = xyz->x(str_ind) - xyz->x(b_ind);
	vec[2] = xyz->y(str_ind) - xyz->y(b_ind);
	vec[3] = xyz->z(str_ind) - xyz->z(b_ind);
	vec[0] = sqrt( (vec[1]*vec[1]) + (vec[2]*vec[2]) + (vec[3]*vec[3]) );
	double rat_bb = b_len / vec[0]; //ratio of bond lengths for old and new
	//scale vector by this ratio
	vec[1] = vec[1]*rat_bb;
	vec[2] = vec[2]*rat_bb;
	vec[3] = vec[3]*rat_bb;
	//add this to the bonding atoms position to get the capping atom
	wrtnum++;
	ext_crd->sett(wrtnum, a_type);
	double n_pos = xyz->x(b_ind) + vec[1];
	ext_crd->setx(wrtnum, n_pos);
	n_pos = xyz->y(b_ind) + vec[2];
	ext_crd->sety(wrtnum, n_pos);
	n_pos = xyz->z(b_ind) + vec[3];
	ext_crd->setz(wrtnum, n_pos);
	
	
	int offset = num_init - e_ind;
	if(e_ind < num_init){ //must shift atom indecies if some atoms at end are not extracted
		for(int jj = 1; jj <= offset; jj++){
			wrtnum++; 
			int pre_ind = e_ind + jj ;
			ext_crd->sett(wrtnum, xyz->type(pre_ind) ) ;
			ext_crd->setx(wrtnum, xyz->x(pre_ind) );
			ext_crd->sety(wrtnum, xyz->y(pre_ind) );
		    ext_crd->setz(wrtnum, xyz->z(pre_ind) );
		}
	}

	ext_crd->setn(wrtnum);
return ext_crd;
} 


void translate_o(XYZ_CRD * xyz, int a_ind){
	// translates Origin to coinside with specified atom index.
	double dx = xyz->x(a_ind);
	double dy = xyz->y(a_ind);
	double dz = xyz->z(a_ind);
	for(int ii=1; ii<= xyz->n(); ii++){
		double nx = xyz->x(ii) - dx;
		xyz->setx(ii,nx);
		double ny = xyz->y(ii) - dy;
		xyz->sety(ii,ny);
		double nz = xyz->z(ii) - dz;
		xyz->setz(ii,nz);
	}
}
	
void rotate_x(XYZ_CRD * xyz, double Theta){
	//rotates all atoms around x axis by angle theta
    for(int ii=1; ii<= xyz->n(); ii++){
		double nY=(cos(Theta)* xyz->y(ii) )-(sin(Theta)* xyz->z(ii) );
		double nZ=(sin(Theta)* xyz->y(ii) )+(cos(Theta)* xyz->z(ii) );
		xyz->sety(ii,nY);
		xyz->setz(ii,nZ);
	}
}
	
void rotate_y(XYZ_CRD * xyz, double Theta){
	//rotates all atoms around y axis by angle theta
    for(int ii=1; ii<= xyz->n(); ii++){
		double nZ=(cos(Theta)* xyz->z(ii) )-(sin(Theta)* xyz->x(ii) );
		double nX=(sin(Theta)* xyz->z(ii) )+(cos(Theta)* xyz->x(ii) );
		xyz->setz(ii,nZ);
		xyz->setx(ii,nX);
	}		
}

void rotate_z(XYZ_CRD * xyz, double Theta){
	//rotates all atoms around y axis by angle theta
	for(int ii=1; ii<= xyz->n(); ii++){
		double nX=(cos(Theta)* xyz->x(ii) )-(sin(Theta)* xyz->y(ii) );
		double nY=(sin(Theta)* xyz->x(ii) )+(cos(Theta)* xyz->y(ii) );
		xyz->setx(ii,nX);
		xyz->sety(ii,nY);
	}		
}

//XYZ_VEC_CRD functions

void translate_xyzvec(double dx, double dy, double dz,XYZ_VEC_CRD * xyz_trans){ //translates index atom to origin
	Log(3) <<"Translate system "<<dx<<" in x, "<<dy<<" in y, "<<dz<<" in z\n";
	for(int ii=1;ii<= xyz_trans->n() ;ii++){
		double xdx = xyz_trans->x(ii) + dx;
		double ydy = xyz_trans->y(ii) + dy;
		double zdz = xyz_trans->z(ii) + dz;
		xyz_trans->setx(ii,xdx);
		xyz_trans->sety(ii,ydy);
		xyz_trans->setz(ii,zdz);
		Log(3) <<"  "<<xyz_trans->type(ii) << " moved to "<<xdx<<" "<<ydy<<" "<<zdz<<"\n";
	}
	
return;
}
	
void translate_o_vec(XYZ_VEC_CRD * xyz, int a_ind){
	// translates Origin to coinside with specified atom index.
	double dx = xyz->x(a_ind);
	double dy = xyz->y(a_ind);
	double dz = xyz->z(a_ind);
	for(int ii=1; ii<= xyz->n(); ii++){
		double nx = xyz->x(ii) - dx;
		xyz->setx(ii,nx);
		double ny = xyz->y(ii) - dy;
		xyz->sety(ii,ny);
		double nz = xyz->z(ii) - dz;
		xyz->setz(ii,nz);
	}
}
	
void rotate_x_vec(XYZ_VEC_CRD * xyz, double Theta){
	//rotates all atoms around x axis by angle theta
    for(int ii=1; ii<= xyz->n(); ii++){
		double nY=(cos(Theta)* xyz->y(ii) )-(sin(Theta)* xyz->z(ii) );
		double nZ=(sin(Theta)* xyz->y(ii) )+(cos(Theta)* xyz->z(ii) );
		xyz->sety(ii,nY);
		xyz->setz(ii,nZ);
	}
}
	
void rotate_y_vec(XYZ_VEC_CRD * xyz, double Theta){
	//rotates all atoms around y axis by angle theta
    for(int ii=1; ii<= xyz->n(); ii++){
		double nZ=(cos(Theta)* xyz->z(ii) )-(sin(Theta)* xyz->x(ii) );
		double nX=(sin(Theta)* xyz->z(ii) )+(cos(Theta)* xyz->x(ii) );
		xyz->setz(ii,nZ);
		xyz->setx(ii,nX);
	}		
}

void rotate_z_vec(XYZ_VEC_CRD * xyz, double Theta){
	//rotates all atoms around y axis by angle theta
	for(int ii=1; ii<= xyz->n(); ii++){
		double nX=(cos(Theta)* xyz->x(ii) )-(sin(Theta)* xyz->y(ii) );
		double nY=(sin(Theta)* xyz->x(ii) )+(cos(Theta)* xyz->y(ii) );
		xyz->setx(ii,nX);
		xyz->sety(ii,nY);
	}		
}


// Torder1 class functions


Torder1::Torder1(){
	x=0; y=0; z=0;
}
Torder1::Torder1(double x1, double y1, double z1){
	x=x1; y=y1; z=z1;
}

double Torder1::mag(){
	double M = x*x;
	M += y*y;
	M += z*z;
return sqrt(M);
}

bool Torder1::isnull(){ //check for null vector
	if(fabs(mag()) < 1e-15 ) return true;
return false;
}

void Torder1::norm(){
	double M = mag();
	if(M==0) return;
	x = x/M;
	y = y/M;
	z = z/M;
return;
}

void Torder1::scale(double a){
	x = x * a;
	y = y * a;
	z = z * a;
return;
}

Torder1& Torder1::operator=(const Torder1 &rhs){
	if(this != &rhs){
		x = rhs.x;
		y = rhs.y;
		z = rhs.z;		
	}
return *this;
}

Torder1& Torder1::operator+=(const Torder1 &rhs){
	x = x + rhs.x;
	y = y + rhs.y;
	z = z + rhs.z;
return *this;
}

Torder1& Torder1::operator-=(const Torder1 &rhs){
	x = x - rhs.x;
	y = y - rhs.y;
	z = z - rhs.z;
return *this;
}

Torder1& Torder1::operator%=(double a){  //scale vector
	x = x * a;
	y = y * a;
	z = z * a;
return *this;
}


Torder1& Torder1::operator|=(const Torder1 &rhs){// cross product
		//cout << "("<<x<<","<<y<<","<<z<<")x("<<rhs.x<<","<<rhs.y<<","<<rhs.z<<")";
	double bx,by,bz;
	bx = (y * rhs.z) - (z * rhs.y);
	by =-(x * rhs.z) + (z * rhs.x);
	bz = (x * rhs.y) - (y * rhs.x);
	x=bx; y=by; z=bz;
		//cout <<" = ("<<x<<","<<y<<","<<z<<")\n";
return *this;
}

// 2 part operators

Torder1 Torder1::proj(const Torder1 &v){ //proj_u(v) = (u dot v)u/(u dot u)
	Torder1 P;
	double uv = dotp(v,*this);
	double uu = (x*x)+(y*y)+(z*z);
	if(uu < 1e-15) return P;
	double a=uv/uu;
	P.x = x*a;
	P.y = y*a;
	P.z = z*a;
	return P;
}

const Torder1 Torder1::operator+(const Torder1 &other) const{
   Torder1 result = *this;  //make copy
   result += other;         //add other to copy
return result;
}

const Torder1 Torder1::operator-(const Torder1 &other) const{
   Torder1 result = *this;  //make copy
   result -= other;         //sub other to copy
return result;
}

const Torder1 Torder1::operator|(const Torder1 &other) const{
   Torder1 result = *this;  //make copy
   result |= other;         //cross other to copy
return result;
}

ostream& operator<<(ostream& os, const Torder1& P){
	os << "("<< P.x <<","<< P.y <<","<< P.z <<")";
return os;
}

/* old dot product
double Torder1::operator*=(const Torder1 &rhs){// dot product
	//Torder1 prod;
	cout <<*this<<" * "<<rhs<<"";
	double dot = x * rhs.x;
	dot += y * rhs.y;
	dot += z * rhs.z;
	cout <<" = "<<dot<<"\n\n";
return dot;
}
const Torder1 Torder1::operator*(const Torder1 &other) const{
   Torder1 result = *this;  //make copy
   result *= other;         //dot other to copy
return result;
}

*/

Torder1 operator%(double a,const Torder1 &rhs){
	Torder1 mult(rhs.x,rhs.y,rhs.z);
	mult %= a;
return mult;
}
Torder1 operator%(const Torder1 &lhs,double a){
	Torder1 mult(lhs.x,lhs.y,lhs.z);
	mult %= a;
return mult;
}

double dotp(const Torder1& v1, const Torder1& v2 ){
	//cout <<v1<<" * "<<v2<<"";
	double dot = v1.x * v2.x;
	dot += v1.y * v2.y;
	dot += v1.z * v2.z;
	//cout <<" = "<<dot<<"\n\n";
return dot;
}

//other functions
Torder1 getAP(XYZ_VEC_CRD * xyz, int a_ind){
	Torder1 P;
	P.x = xyz->x(a_ind);
	P.y = xyz->y(a_ind);
	P.z = xyz->z(a_ind);
return P;
}

