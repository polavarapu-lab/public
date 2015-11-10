#ifndef XYZ_MANIP_H_
#define XYZ_MANIP_H_


XYZ_CRD * strip_mol(XYZ_CRD * xyz, int srt_ind, int e_ind, int b_ind, double b_len, string a_type);
void translate_o(XYZ_CRD * xyz, int a_ind);   // translates Origin to coinside with specified atom index.
void rotate_x(XYZ_CRD * xyz, double theta);   //rotates all atoms around x axis by angle theta
void rotate_y(XYZ_CRD * xyz, double theta);   //rotates all atoms around y axis by angle theta
void rotate_z(XYZ_CRD * xyz, double theta);   //rotates all atoms around z axis by angle theta

//xyzvec types
void translate_xyzvec(double dx, double dy, double dz,XYZ_VEC_CRD * xyz_trans);
void translate_o_vec(XYZ_VEC_CRD * xyz, int a_ind);
void rotate_x_vec(XYZ_VEC_CRD * xyz, double Theta);
void rotate_y_vec(XYZ_VEC_CRD * xyz, double Theta);
void rotate_z_vec(XYZ_VEC_CRD * xyz, double Theta);

//coordinate class for doing manipulations
class Torder1 { //tensor rank 1, 3 components
	public:
	double x;
	double y;
	double z;
	Torder1();
	Torder1(double x1,double y1, double z1);
	double mag(); //returns the magnitude of the vector
	bool isnull(); //returns true if magnitude is less than 1e-15
	//operator functions on self
	void norm(); //norm vector
	void scale(double a);
	Torder1& operator=(const Torder1 &rhs);            //copy
	Torder1& operator+=(const Torder1 &rhs);           //add to
	Torder1& operator-=(const Torder1 &rhs);           //subtract from
    Torder1& operator%=(double a);                     //scale vector by a
		//Torder1& operator*=(const Torder1 &rhs);       //dot with
	Torder1& operator|=(const Torder1 &rhs);           //cross with
	
	
	//2 part operators
	Torder1 proj(const Torder1 &v); // projection of v in direction of u
	const Torder1 operator+(const Torder1 &rhs) const;   //add      
	const Torder1 operator-(const Torder1 &rhs) const;   //subtract
		//double operator*(const Torder1 &rhs) const;   //dot
	const Torder1 operator|(const Torder1 &rhs) const;   //cross
	
	//output operator
	friend ostream& operator<<(ostream& os, const Torder1& P);
	
//see  http://courses.cms.caltech.edu/cs11/material/cpp/donnie/cpp-ops.html  for operator definitions
} ;
//defined outside of class to avoid implied *this
	Torder1 operator%(double a,const Torder1 &rhs);   //scale (a*v)
	Torder1 operator%(const Torder1 &lhs,double a);   //scale (v*a)
//other functions
double dotp(const Torder1& v1, const Torder1& v2 );

// functions with other classes
Torder1 getAP(XYZ_VEC_CRD * xyz, int a_ind); //get atom index 

#endif


