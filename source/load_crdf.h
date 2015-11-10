

#ifndef LOAD_CRDF_H_
#define LOAD_CRDF_H_

#include <vector>

class XYZ_CRD {
	private:
		int num_l; //number of crd stored
		double psx[6000];
		double psy[6000];
		double psz[6000];
		string atom_type[6000];
		string file_comment;
	public:
		XYZ_CRD();
		double x(int BIN); //get x crd
		void setx(int BIN, double VAL); //set x crd
		double y(int BIN);
		void sety(int BIN, double VAL);
		double z(int BIN);
		void setz(int BIN, double VAL);
		int n(); //get number stored
		void setn(int VAL);
		string type(int BIN);
		void sett(int BIN,string VAL);
		string comment(); //get the comment line
		void set_comment(string comm1); 

} ;

class XYZ_VEC_CRD {
	private:
		int num_l; //number of crd stored
		int num_MM;
		std::vector<double> psx;
		std::vector<double> psy;
		std::vector<double> psz;
		std::vector<double> chr;
		std::vector<double> pol;
		std::vector<string> atom_type;
		std::vector<string> res;
		string file_comment;
		string desc; //description
		double gendat[9]; //general data fields
	public:
		XYZ_VEC_CRD();
		XYZ_VEC_CRD(string descript);
		double x(int BIN); //get x crd
		void setx(int BIN, double VAL); //set x crd
		double y(int BIN);
		void sety(int BIN, double VAL);
		double z(int BIN);
		void setz(int BIN, double VAL);
		double c(int BIN);
		void setc(int BIN, double VAL);	
		double p(int BIN);
		void setp(int BIN, double VAL);	
		double dat(int BIN);
		void setdat(int BIN, double VAL);
		
		int n(); //get number of xyz crd stored
		void setn(int VAL);
		int MM(); //get number of zyxchr stored
		void setMM(int VAL);
		string gres(int BIN);
		void setr(int BIN,string VAL);
		string type(int BIN);
		void sett(int BIN,string VAL);
		string comment(); //get the comment line
		void set_comment(string comm1); 
		string getdesc(); //retrn description
		void set_desc(string DES1); 
		
		// simple analysis functions
		bool is_Hydrogen(int BIN);      // returns true if atom type is hydrogen -- H or h
		bool is_Oxygen(int BIN);        // returns true if atom type is oxygen   -- O or o
		bool is_hydroxyl(int BIN,double mdist); // detects if atom at BIN is member of a hydroxyl group. 		
		
		//loading functions
		bool add_pos(string TP, double X, double Y, double Z); //add atom Do not use setx, sety, setz to add new crd. 
		bool add_pcp(string TP, double X, double Y, double Z, double C, double P, string RE); 
		bool copy_pcp( XYZ_VEC_CRD * xyzp, int idx); //adds xyzchr new atom from xyzp crds at index idx
		   //check to see if there are charge and pol info before copying (both xyz and xyzchr)
		bool copy_dual( XYZ_VEC_CRD * xyzp, int idx); //adds new atom from xyzp crds at index idx
		bool copy_info( XYZ_VEC_CRD * xyzp); // copies all comment, gendat, and desc from xyzp to new set. 
		
		// other misc functions
		void print();                    // print coordinates to screen
		
		//clearing functions
		bool rem_pos(int idx); //removes indexed atom from xyz_vec_crd
		void cclear();
		
} ;

class RES_XYZ {
	private:
		int num_l; //number of crd stored
		double psx[6000];
		double psy[6000];
		double psz[6000];
		double property1[6000];
		double property2[6000];
		string atom_type[6000];
		string atom_res[6000];
		string file_comment;
	public:
		RES_XYZ();
		double x(int BIN); //get x crd
		void setx(int BIN, double VAL); //set x crd
		double y(int BIN);
		void sety(int BIN, double VAL);
		double z(int BIN);
		void setz(int BIN, double VAL);
		int n(); //get number stored
		void setn(int VAL);
		string type(int BIN);
		void sett(int BIN,string VAL);
		string res(int BIN);
		void setres(int BIN, string res);
		string comment(); //get the comment line
		void set_comment(string comm1); 
		void setprop1(int BIN, double p1);
		double prop1(int BIN);
		void setprop2(int BIN, double p2);
		double prop2(int BIN);


} ;

XYZ_CRD * loadcrd_mol(int MAX, string filename);
XYZ_CRD * loadcrd_pdb(int MAX, string filename,bool p2_flag);
XYZ_CRD * loadcrd_xyz(int MAX, string filename); 
int printcrd_xyz(XYZ_CRD * xyzp_crd, string filename1, string comment1); //print xyz coordinates


string AssignType(int AtomN); //gives atomic symbol from atomic number
double AssignMass(string type);
int type2atomN(string type);  //gives atomic number from symbol

// string manip
double StoF(string CONV);
inline std::string trimstr(std::string& str);

// other

XYZ_VEC_CRD * loadcrd_pdbvec(int MAX, string filename,bool p2_flag);
XYZ_VEC_CRD * loadcrd_xyzvec(int MAX, string filename);
XYZ_VEC_CRD * loadcrd_xyzchr_vec(int MAX, string filename);
XYZ_VEC_CRD * loadcrd_molvec(int MAX, string filename);
XYZ_VEC_CRD * loadcrd_amberRST(int Max, string topFN, string rstFN);

int printcrd_xyzvec(XYZ_VEC_CRD * xyzp_crd, string filename1, string comment1); //print crd to file
int printcrd_xyzvec_precise(XYZ_VEC_CRD * xyzp_crd, string filename1, string comment1); // print crd with high precision
int printcrd_xyzchr_vec(XYZ_VEC_CRD * xyzp_crd, string filename1, string comment1); //print xyzchr crd to file
int printcrd_xyzvec_range(XYZ_VEC_CRD * xyzp_crd,int str,int fin, string filename1, string comment1); //print crd in range

/*exceptions
99 = file load.
88 = xyzchr number of xyz does not match

*xyzchr file format
line1=number of atoms (natom)
line2=comment
lines 3 to natom+2 = xyz crd with charges and polarizabilities 
	/Format= Atomic_Symbol  X_crd  Y_crd  Z_crd  Charge  Polarizability  Residue_Label
Line natom+3=Blank
Line natom+4="--DATA--"
	3 lines of atomic boundary conditions (xyz)
	3 lines of atomic boundary conditions (angles)
	3 lines (not yet used)


*/
#endif
