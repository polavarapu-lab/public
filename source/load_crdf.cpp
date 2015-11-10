#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <math.h>
#include <vector>
#include <iomanip>
#include <exception>
using namespace std;

#include "load_crdf.h"
#include "class_log.h" //use of Log() function
/*
class XYZ_CRD {
	private:
		int num_l; //number of crd stored
		double psx[6000];
		double psy[6000];
		double psz[6000];
		string atom_type[6000];
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

} ;*/
//XYZ_CRD class functions

XYZ_CRD::XYZ_CRD(){
	num_l=0;
	for(int ii=0; ii<6000; ii++){
		psx[ii]=0;
		psy[ii]=0;
		psx[ii]=0;
		atom_type[ii]="";
	}
	file_comment="";
}
double XYZ_CRD::x(int BIN){
	return psx[BIN];
}
double XYZ_CRD::y(int BIN){
	return psy[BIN];
}
double XYZ_CRD::z(int BIN){
	return psz[BIN];
}
int XYZ_CRD::n(){
	return num_l;
}
string XYZ_CRD::type(int BIN){
	return atom_type[BIN];
}

void XYZ_CRD::setx(int BIN,double VAL){
	psx[BIN]=VAL;
return;
}
void XYZ_CRD::sety(int BIN,double VAL){
	psy[BIN]=VAL;
return;
}
void XYZ_CRD::setz(int BIN,double VAL){
	psz[BIN]=VAL;
return;
}
void XYZ_CRD::setn(int VAL){
	num_l=VAL;
return;
}
void XYZ_CRD::sett(int BIN,string VAL){
	atom_type[BIN]=VAL;
return;
}
string XYZ_CRD::comment(){
return file_comment;
}
void XYZ_CRD::set_comment(string comm1){
	file_comment = comm1;
return;
}

RES_XYZ::RES_XYZ(){
	num_l=0;
	for(int ii=0; ii<6000; ii++){
		psx[ii]=0;
		psy[ii]=0;
		psx[ii]=0;
		property1[ii]=0;
		property2[ii]=0;
		atom_type[ii]="";
		atom_res[ii]="";
	}
	file_comment="";
}
double RES_XYZ::x(int BIN){
	return psx[BIN];
}
double RES_XYZ::y(int BIN){
	return psy[BIN];
}
double RES_XYZ::z(int BIN){
	return psz[BIN];
}
int RES_XYZ::n(){
	return num_l;
}
string RES_XYZ::type(int BIN){
	return atom_type[BIN];
}

void RES_XYZ::setx(int BIN,double VAL){
	psx[BIN]=VAL;
return;
}
void RES_XYZ::sety(int BIN,double VAL){
	psy[BIN]=VAL;
return;
}
void RES_XYZ::setz(int BIN,double VAL){
	psz[BIN]=VAL;
return;
}
void RES_XYZ::setn(int VAL){
	num_l=VAL;
return;
}
void RES_XYZ::sett(int BIN,string VAL){
	atom_type[BIN]=VAL;
return;
}
string RES_XYZ::comment(){
return file_comment;
}
void RES_XYZ::set_comment(string comm1){
	file_comment = comm1;
return;
}
void RES_XYZ::setprop1(int BIN, double p1){
	property1[BIN] = p1;
return;
}
double RES_XYZ::prop1(int BIN){
	return property1[BIN];
}
void RES_XYZ::setprop2(int BIN, double p2){
	property2[BIN] = p2;
return;
}
double RES_XYZ::prop2(int BIN){
	return property2[BIN];
}
void RES_XYZ::setres(int BIN,string res){
	atom_res[BIN]=res;
return;
}
string RES_XYZ::res(int BIN){
return atom_res[BIN];
}

double StoF(string CONV){

  double val ;
  stringstream ss (stringstream::in | stringstream::out);
  ss << CONV;
  ss >> val ;
  //cout << "float" << val <<endl;
  return val;
}

XYZ_CRD * loadcrd_mol(int MAX, string filename){
	XYZ_CRD * xyzmol = new XYZ_CRD;
	ifstream infile;
    infile.open (filename.c_str());
	if(infile.fail()==true){ 
		cerr << "Input File "<<filename<<" Not Found" << endl; 
		infile.close();
		throw 99;
		return xyzmol; 
	}
	Log(1) <<" input file .mol selected \n"; 
	string STRING,LINE1;
	getline(infile,STRING); getline(infile,STRING); //get 4 comment lines for mol fies
	getline(infile,STRING); getline(infile,STRING);	
		if(MAX > 5999){
		MAX = 5999; cerr << "max read set at 5999"<<endl;
	}
	for(int ii=1;ii<=MAX;ii++){
		getline(infile,LINE1);
		if(LINE1.length() < 40){ //end of file, or error in file or line
			infile.close();
			return xyzmol;
		}
        xyzmol->sett(ii,LINE1.substr(31,2) ); //get atomtype
        xyzmol->setx(ii, StoF(LINE1.substr(0,11)) );   
        xyzmol->sety(ii, StoF(LINE1.substr(12,9)) );
        xyzmol->setz(ii, StoF(LINE1.substr(22,9)) );
		xyzmol->setn(ii);
		// cout << "atom added, type " << xyzmol->type(ii) << "  "<< xyzmol->x(ii) << "  "<< xyzmol->y(ii)<< "  "<< xyzmol->z(ii)<< "  NUM "<< xyzmol->n() << endl;
		if(!infile.eof()==0){  
			cerr << "possible error in coordinate file" << endl;
			infile.close();
			return xyzmol; 
		}
	}
	infile.close();
	if(xyzmol->n() == MAX){
		cerr <<"max number of atoms reached."<< endl;
	}
return xyzmol;
}

XYZ_CRD * loadcrd_pdb(int MAX, string filename,bool p2_flag){
	XYZ_CRD * xyzpdb = new XYZ_CRD;
	ifstream infile;
    infile.open (filename.c_str());
	if(infile.fail()==true){ 
		cerr << "Input File "<<filename<<" Not Found" << endl; 
		infile.close();
		throw 99;
		return xyzpdb; 
	}
	Log(1) <<" input file .pdb selected \n";  
	//variable number of comment lines in pdb
	string LINE1;
	if(MAX > 5999){
		MAX = 5999; cerr << "max read set at 5999"<<endl;
	}
	for(int ii=1;ii<=MAX;){
		getline(infile,LINE1);
		if(LINE1.substr(0,3)=="END"){
			infile.close();
			return xyzpdb; 
		}
		if(!infile.eof()==0){
			infile.close();
			return xyzpdb;
		}
		if(LINE1.substr(0,4)!="ATOM") continue;
		if(LINE1.length() < 60){ //error in file or line
			cout <<"stopped after short line"<<endl;
			infile.close();
			return xyzpdb; 
		} 
        string a2 = LINE1.substr(13,1);
		if(p2_flag==false){
			string a3 = LINE1.substr(14,1);
			if(a3=="1"){}      else if(a3=="2"){}
			else if(a3=="3"){} else if(a3=="4"){}
			else if(a3=="5"){} else if(a3=="6"){}
			else if(a3=="7"){} else if(a3=="8"){}
			else if(a3=="9"){} else if(a3=="0"){}
			else{ a2 += a3;  } 
		}
			//a2.Trim(); cout << "Type detected "<<a2<<endl;
        xyzpdb->sett(ii, a2);
		//cout << "x = '"<<LINE1.substr(30,8)<<"', ";
		//cout << "x = '"<<LINE1.substr(38,8)<<"', ";
		//cout << "x = '"<<LINE1.substr(46,8)<<"'\n";
        xyzpdb->setx(ii, StoF(LINE1.substr(30,8)) );   
        xyzpdb->sety(ii, StoF(LINE1.substr(38,8)) );
        xyzpdb->setz(ii, StoF(LINE1.substr(46,8)) );
		xyzpdb->setn(ii);
		string check2=LINE1.substr(30,30);
		double comp[4];
		std::istringstream iss(check2);
		iss >> comp[1]; 
		iss >> comp[2];
		iss >> comp[3];
		if((xyzpdb->x(ii) - comp[1]) > 0.01){
			cerr <<" Warning check x pos atom "<<ii<<" got "<<comp[1]<<" and "<<xyzpdb->x(ii)<<" from string "<<LINE1<<endl;
		}
		if((xyzpdb->y(ii) - comp[2]) > 0.01){
			cerr <<" Warning check y pos atom "<<ii<<" got "<<comp[2]<<" and "<<xyzpdb->y(ii)<<" from string "<<LINE1<<endl;
		}

		if((xyzpdb->z(ii) - comp[3]) > 0.01){
			cerr <<" Warning check z pos atom "<<ii<<" got "<<comp[3]<<" and "<<xyzpdb->z(ii)<<" from string "<<LINE1<<endl;
		}
		ii++; //atom index advances from successful loading
		//cout << "atom added, type " << atom[CNT1].type << " pol "<< atom[CNT1].bondPOL << endl;
		//cout <<","<<LINE1.substr(31,7)<<","<<LINE1.substr(39,7)<<","<<LINE1.substr(47,7)<<endl;
	}
	infile.close();
	if(xyzpdb->n() == MAX){
		cerr <<"max number of atoms reached."<< endl;
	}
return xyzpdb;
}

XYZ_CRD * loadcrd_xyz(int MAX, string filename){
	XYZ_CRD * xyz = new XYZ_CRD;
	ifstream infile;
    infile.open (filename.c_str());
	if(infile.fail()==true){ 
		cerr << "Input File "<<filename<<" Not Found" << endl; 
		infile.close();
		throw 99;
		return xyz; 
	}
	Log(1) <<" input file .xyz selected \n";  
	//read comment lines in xyz
	string LINE1;
	getline(infile,LINE1);
	int num_af = (int)StoF(LINE1); //read number of atoms in file from line1
	getline(infile,LINE1);
	xyz->set_comment(LINE1);
	if(MAX > 5999){
		MAX = 5999; cerr << "max read set at 5999"<<endl;
	}
	for(int ii=1;ii<=MAX;ii++){
		getline(infile,LINE1);
		string FRAG;
		stringstream ss(LINE1);
		vector<string> tokens;
		int cc=1;
		while( ss >> FRAG){
			tokens.push_back(FRAG);
			if(cc==1) xyz->sett(ii,FRAG) ; //atomtype is 1st 
			else if(cc==2) xyz->setx(ii,StoF(FRAG) ); //x coord 
			else if(cc==3) xyz->sety(ii,StoF(FRAG) ); //y coord
			else if(cc==4) xyz->setz(ii,StoF(FRAG) ); //z coord
			cc++;
		}
		if(cc<5){ 
			infile.close();
			if(xyz->n() != num_af){
				cerr <<"Number of atoms found "<<xyz->n()<<" is different from number specified in file "<<num_af<<endl;
			}
			return xyz; //error in file or line
		}
		else{
			xyz->setn(ii);
		}
	}
infile.close();
	if(xyz->n() == MAX){
		cerr <<"max number of atoms reached."<< endl;
	}
return xyz;
}
int printcrd_xyz(XYZ_CRD * xyzp_crd, string filename1, string comment1){
	ofstream outfile2 ( filename1.c_str() );
	int CNT6 = xyzp_crd->n();
	outfile2 << CNT6 << endl;
	outfile2 << comment1 << endl;
	outfile2 << fixed <<showpoint;
	//outfile2 << setprecision(5);
	char XP[14],YP[14],ZP[14];
	for( int xx=1; xx<= CNT6; xx++){ // print xyz coor
			sprintf(XP,"%12.6f",xyzp_crd->x(xx) );
			sprintf(YP,"%12.6f",xyzp_crd->y(xx) );
			sprintf(ZP,"%12.6f",xyzp_crd->z(xx) );
			//outfile2  << xyzp_crd->type(xx) << "\t"  << xyzp_crd->x(xx) << "\t" << xyzp_crd->y(xx) << "\t" << xyzp_crd->z(xx) << endl;
			outfile2 << xyzp_crd->type(xx) << "\t"  << XP << "\t"  << YP << "\t"  << ZP <<endl;
	}
	outfile2.close();
	Log() <<"Wrote "<<CNT6<<" to "<<filename1<<"\n";
return 0;
}
/*

bool is_number(char c_a){ //check to see if a character is a number 0-9
	if(c_a=='1') return true;
	if(c_a=='2') return true;
	if(c_a=='3') return true;
	if(c_a=='4') return true;
	if(c_a=='5') return true;
	if(c_a=='6') return true;
	if(c_a=='7') return true;
	if(c_a=='8') return true;
	if(c_a=='9') return true;
	if(c_a=='0') return true;
return false;
}

*/


//XYZ_VEC_CRD class functions

XYZ_VEC_CRD::XYZ_VEC_CRD(){
	num_l=0;
	num_MM = 0;
	file_comment="";
	psx.reserve(100000);
	psy.reserve(100000);
	psz.reserve(100000);
	chr.reserve(100000);
	pol.reserve(100000);
	atom_type.reserve(100000);
    res.reserve(100000);
	desc="";
}

XYZ_VEC_CRD::XYZ_VEC_CRD(string descript){
	num_l=0;
	num_MM = 0;
	file_comment="";
	psx.reserve(100000);
	psy.reserve(100000);
	psz.reserve(100000);
	chr.reserve(100000);
	pol.reserve(100000);
	atom_type.reserve(100000);
    res.reserve(100000);
	desc=descript;
	for(int i=0;i<9;i++) gendat[i]=0;
}

double XYZ_VEC_CRD::x(int BIN){
   try{	return psx.at(BIN-1); }
   catch(std::exception& o){ std::cerr<<o.what()<<" "<<desc<<" read psx i="<<BIN<<std::endl; }
	cerr << "error bad read X number "<<BIN<<"\n";
return 0;
}
double XYZ_VEC_CRD::y(int BIN){
   try{	return psy.at(BIN-1); }
   catch(std::exception& o){ std::cerr<<o.what()<<" "<<desc<<" read psy i="<<BIN<<std::endl; }
	cerr << "error bad read Y number "<<BIN<<"\n";
return 0;
}
double XYZ_VEC_CRD::z(int BIN){
   try{	return psz.at(BIN-1); }
   catch(std::exception& o){ std::cerr<<o.what()<<" "<<desc<<" read psz i="<<BIN<<std::endl; }
	cerr << "error bad read Z number "<<BIN<<"\n";
return 0;
}
double XYZ_VEC_CRD::c(int BIN){
   try{	return chr.at(BIN-1); }
   catch(std::exception& o){ std::cerr<<o.what()<<" "<<desc<<" read chr i="<<BIN<<std::endl; }
	cerr << "error bad read charge number "<<BIN<<"\n";
return 0;
}
double XYZ_VEC_CRD::p(int BIN){
   try{	return pol.at(BIN-1); }
   catch(std::exception& o){ std::cerr<<o.what()<<" "<<desc<<" read pol i="<<BIN<<std::endl; }
	cerr << "error bad read charge number "<<BIN<<"\n";
return 0;
}
double XYZ_VEC_CRD::dat(int BIN){
	if(BIN > 9 || BIN < 1) cerr << "error bad read from general data array "<<BIN<<". Only 1-9 possible.\n";
return gendat[BIN-1];
}
int XYZ_VEC_CRD::n(){
	return num_l;
}
int XYZ_VEC_CRD::MM(){
	return num_MM;
}
string XYZ_VEC_CRD::type(int BIN){
   try{	return atom_type.at(BIN-1); }
   catch(std::exception& o){ std::cerr<<o.what()<<" "<<desc<<std::endl; }
   cerr <<" "<<desc<< " error bad read of atom type\n";
return ""; 
}
string XYZ_VEC_CRD::gres(int BIN){
   try{	return res.at(BIN-1); }
   catch(std::exception& o){ std::cerr<<o.what()<<" "<<desc<<std::endl; }
   cerr <<" "<<desc<<" error bad read of ResID\n";
return ""; 
}

void XYZ_VEC_CRD::setx(int BIN, double VAL){
   try{	 psx.at(BIN-1) = VAL; }
   catch(std::exception& o){ std::cerr<<o.what()<<" "<<desc<<" write psx i="<<BIN<<std::endl; }
return;
}
void XYZ_VEC_CRD::sety(int BIN, double VAL){
   try{	 psy.at(BIN-1) = VAL; }
   catch(std::exception& o){ std::cerr<<o.what()<<" "<<desc<<" write psy i="<<BIN<<std::endl; }
return;
}
void XYZ_VEC_CRD::setz(int BIN, double VAL){
   try{	 psz.at(BIN-1) = VAL; }
   catch(std::exception& o){ std::cerr<<o.what()<<" "<<desc<<" write psz i="<<BIN<<std::endl; }
return;
}
void XYZ_VEC_CRD::setc(int BIN, double VAL){
   try{	 chr.at(BIN-1) = VAL; }
   catch(std::exception& o){ std::cerr<<o.what()<<" "<<desc<<" write chr i="<<BIN<<std::endl; }
return;
}
void XYZ_VEC_CRD::setp(int BIN, double VAL){
   try{	 pol.at(BIN-1) = VAL; }
   catch(std::exception& o){ std::cerr<<o.what()<<" "<<desc<<" write pol i="<<BIN<<std::endl; }
return;
}
void XYZ_VEC_CRD::setdat(int BIN,double VAL){
	if(BIN > 9 || BIN < 1) cerr << "error bad set for general data array "<<BIN<<". Only 1-9 possible.\n";
	gendat[BIN-1] = VAL;
return;
}

void XYZ_VEC_CRD::setn(int VAL){
	num_l=VAL;
return;
}
void XYZ_VEC_CRD::setMM(int VAL){
	num_MM=VAL;
return;
}
void XYZ_VEC_CRD::sett(int BIN, string VAL){
   try{	 atom_type.at(BIN-1) = VAL; }
   catch(std::exception& o){ std::cerr<<o.what()<<" "<<desc<<std::endl; }
return;
}
void XYZ_VEC_CRD::setr(int BIN, string VAL){
   try{	 res.at(BIN-1) = VAL; }
   catch(std::exception& o){ std::cerr<<o.what()<<" "<<desc<<std::endl; }
return;
}
string XYZ_VEC_CRD::comment(){
return file_comment;
}

void XYZ_VEC_CRD::set_comment(string comm1){
	file_comment = comm1;
return;
}

string XYZ_VEC_CRD::getdesc(){
return desc;
}

void XYZ_VEC_CRD::set_desc(string DES1){
	desc = DES1;
return;
}

bool XYZ_VEC_CRD::is_Hydrogen(int BIN){    // returns true if atom type is hydrogen -- H or h
	string chk_type; 
   try{	 chk_type = atom_type.at(BIN-1);  }
   catch(std::exception& o){ std::cerr<<o.what()<<" "<<desc<<std::endl; }
   trimstr(chk_type);
   if( chk_type.length() > 1 ) return false;
   if( chk_type == "H" ) return true;
   if( chk_type =="h" ) return true;
return false;
}

bool XYZ_VEC_CRD::is_Oxygen(int BIN){    // returns true if atom type is oxygen   -- O or o
	string chk_type; 
   try{	 chk_type = atom_type.at(BIN-1);  }
   catch(std::exception& o){ std::cerr<<o.what()<<" "<<desc<<std::endl; }
   trimstr(chk_type);
   if( chk_type.length() > 1 ) return false;
   if( chk_type == "O" ) return true;
   if( chk_type =="o" ) return true;
return false;
}

bool XYZ_VEC_CRD::is_hydroxyl(int BIN,double mdist){  // detects if atom at BIN is member of a hydroxyl group. 
 // must be H or O and must have other atom within mdist Angstroms
 double smallest_d=9e99; // smallest d for bond partner
 if( is_Hydrogen(BIN) ){
    for( int ii=1; ii<= num_l; ii++){
	   if(ii==BIN) continue;
	   // calc dist for all to check for small distances
	   double x = psx.at(BIN-1) - psx.at(ii-1);
	   double y = psy.at(BIN-1) - psy.at(ii-1);
	   double z = psz.at(BIN-1) - psz.at(ii-1);
	   double d = sqrt( x*x + y*y + z*z );
	   if( d < 1e-4 ){
	      cerr << "warning small distance "<<d<<" found. Check coordinates.\n bool XYZ_VEC_CRD::is_hydroxyl(int BIN,double mdist) failed\n";
		  return false;
	   }
	   if( is_Oxygen(ii) ){ // 
	      if( d < smallest_d ) smallest_d = d;
	   }
    }
	if ( smallest_d < mdist ) return true;
 }
 else if ( is_Oxygen(BIN) ){
    for( int ii=1; ii<= num_l; ii++){
	   if(ii==BIN) continue;
	   // calc dist for all to check for small distances
	   double x = psx.at(BIN-1) - psx.at(ii-1);
	   double y = psy.at(BIN-1) - psy.at(ii-1);
	   double z = psz.at(BIN-1) - psz.at(ii-1);
	   double d = sqrt( x*x + y*y + z*z );
	   if( d < 1e-4 ){
	      cerr << "warning small distance "<<d<<" found. Check coordinates.\n bool XYZ_VEC_CRD::is_hydroxyl(int BIN,double mdist) failed\n";
		  return false;
	   }
	   if( is_Hydrogen(ii) ){ // 
	      if( d < smallest_d ) smallest_d = d;
	   }
    }
	if ( smallest_d < mdist ) return true; 
 }
return false;
}

bool XYZ_VEC_CRD::add_pos(string TP, double X, double Y, double Z){
	//check to make sure all vectors have the same number of elements
	int cn = psx.size();
	if(cn != num_l){ 
		cerr<<"Number of atoms stored ("<<psy.size()<<") not equal to X elements "<<psx.size()<<"\n";
		return false;
	}
	if(cn != (int)psy.size()){ 
		cerr<<"Number of Y ("<<psy.size()<<") crd not equal to X "<<psx.size()<<"\n";
		return false;
	}
	if(cn != (int)psz.size()){ 
		cerr<<"Number of Z ("<<psz.size()<<") crd not equal to X "<<psx.size()<<"\n";
		return false;
	}
	if(cn != (int)atom_type.size()){ 
		cerr<<"Number of Labels "<<atom_type.size()<<" not equal to X "<<psx.size()<<"\n";
		return false;
	}	
	atom_type.push_back(TP);
	psx.push_back(X);
	psy.push_back(Y);
	psz.push_back(Z);
	num_l++;
return true;
}

bool XYZ_VEC_CRD::add_pcp(string TP, double X, double Y, double Z, double C, double P, string RE){
	//check to make sure all vectors have the same number of elements
	int cn = psx.size();
	if(cn != num_l){ 
		cerr<<"Number of atoms stored ("<<psy.size()<<") not equal to X elements "<<psx.size()<<"\n";
		return false;
	}
	if(cn != (int)psy.size()){ 
		cerr<<"Number of Y ("<<psy.size()<<") crd not equal to X "<<psx.size()<<"\n";
		return false;
	}
	if(cn != (int)psz.size()){ 
		cerr<<"Number of Z ("<<psz.size()<<") crd not equal to X "<<psx.size()<<"\n";
		return false;
	}
	if(cn != (int)atom_type.size()){ 
		cerr<<"Number of Labels "<<atom_type.size()<<" not equal to X "<<psx.size()<<"\n";
		return false;
	}	
	if(cn != (int)chr.size()){ 
		cerr<<"Number of charges "<<chr.size()<<" not equal to X "<<psx.size()<<"\n";
		return false;
	}
	if(cn != (int)pol.size()){ 
		cerr<<"Number of polarizabilities "<<pol.size()<<" not equal to X "<<psx.size()<<"\n";
		return false;
	}
	if(cn != (int)res.size()){ 
		cerr<<"Number of ResID "<<pol.size()<<" not equal to X "<<psx.size()<<"\n";
		return false;
	}
	atom_type.push_back(TP);
	psx.push_back(X);
	psy.push_back(Y);
	psz.push_back(Z);
	chr.push_back(C);
	pol.push_back(P);
	res.push_back(RE);
	num_l++;
	num_MM++;
return true;
}

bool XYZ_VEC_CRD::copy_pcp( XYZ_VEC_CRD * xyzp, int idx){  //adds xyzchr new atom from xyzp crds at index idx
	bool res = add_pcp( xyzp->type(idx),xyzp->x(idx),xyzp->y(idx),xyzp->z(idx),xyzp->c(idx),xyzp->p(idx),xyzp->gres(idx) );
return res;
}

bool XYZ_VEC_CRD::copy_dual( XYZ_VEC_CRD * xyzp, int idx){
	bool res;
	int nMM = xyzp->MM();
	if(nMM == 0){ //only xyz crd
		res = add_pos( xyzp->type(idx),xyzp->x(idx),xyzp->y(idx),xyzp->z(idx) );
	}
	else if(nMM == xyzp->n() ){ //same number of xyz and chr and polz
		res = add_pcp( xyzp->type(idx),xyzp->x(idx),xyzp->y(idx),xyzp->z(idx),
			xyzp->c(idx),xyzp->p(idx),xyzp->gres(idx) );
	}
	else if(nMM != xyzp->n() ){
		Log(-1,"XYZ_VEC_CRD::copy_dual( XYZ_VEC_CRD * xyzp, int idx)", 3) <<"Number of crd and charges not equal?\n";
		throw 88;
		return false;
	}
return res;
}

bool XYZ_VEC_CRD::copy_info( XYZ_VEC_CRD * xyzp){
	file_comment = xyzp->comment();
	desc = xyzp->getdesc();
	for(int ii=1;ii<=9; ii++){
		gendat[ii-1] = xyzp->dat(ii);
	}
return true;
}

void XYZ_VEC_CRD::print(){
	char XP[25],YP[25],ZP[25];
	for( int xx=0; xx< num_l; xx++){ // print xyz coor vector counting starts at 0
			sprintf(XP,"%24.18f", psx.at(xx) );
			sprintf(YP,"%24.18f", psy.at(xx) );
			sprintf(ZP,"%24.18f", psz.at(xx) );
			//outfile2  << xyzp_crd->type(xx) << "\t"  << xyzp_crd->x(xx) << "\t" << xyzp_crd->y(xx) << "\t" << xyzp_crd->z(xx) << endl;
			cout << atom_type.at(xx) << " "  << XP << " "  << YP << " "  << ZP <<endl;
	}


}


bool XYZ_VEC_CRD::rem_pos(int idx){ //removes indexed atom from xyz_vec_crd
	if(num_l>0) num_l--;
	psx.erase(psx.begin() + idx);
	psy.erase(psy.begin() + idx);
	psz.erase(psz.begin() + idx);
	atom_type.erase(atom_type.begin() + idx);
	if(num_MM > 0 ){
		chr.erase(chr.begin() + idx);
	    pol.erase(pol.begin() + idx);
	    res.erase(res.begin() + idx);
		num_MM--;
	}
return true;
}
void XYZ_VEC_CRD::cclear(){
	num_l=0;
	num_MM=0;
	atom_type.clear();
	psx.clear();
	psz.clear();
	chr.clear();
	pol.clear();
	res.clear();
	file_comment="";
	desc="";
	for(int i=0;i<9;i++) gendat[i]=0;
return;
}

XYZ_VEC_CRD * loadcrd_pdbvec(int MAX, string filename,bool p2_flag){
	XYZ_VEC_CRD * xyzpdb = new XYZ_VEC_CRD;
	ifstream infile;
    infile.open (filename.c_str());
	if(infile.fail()==true){ 
		cerr << "Input File "<<filename<<" Not Found" << endl; 
		infile.close();
		throw 99;
		return xyzpdb; 
	}
	Log(1) <<" input file .pdb selected \n";  
	//variable number of comment lines in pdb
	string LINE1;
	for(int ii=1;ii<=MAX;){
		getline(infile,LINE1);
		if(LINE1.substr(0,3)=="END"){
			infile.close();
			return xyzpdb; 
		}
		if(!infile.eof()==0){
			infile.close();
			return xyzpdb;
		}
		if(LINE1.substr(0,4)!="ATOM") continue;
		if(LINE1.length() < 60){ //error in file or line
			cout <<"stopped after short line"<<endl;
			infile.close();
			return xyzpdb; 
		} 
        string a2 = LINE1.substr(13,1);
		if(p2_flag==false){
			string a3 = LINE1.substr(14,1);
			if(a3=="1"){}      else if(a3=="2"){}
			else if(a3=="3"){} else if(a3=="4"){}
			else if(a3=="5"){} else if(a3=="6"){}
			else if(a3=="7"){} else if(a3=="8"){}
			else if(a3=="9"){} else if(a3=="0"){}
			else{ a2 += a3;  } 
		}
		xyzpdb->add_pos(a2,StoF(LINE1.substr(30,8)) , StoF(LINE1.substr(38,8)), StoF(LINE1.substr(46,8)) );
		
		string check2=LINE1.substr(30,30);
		double comp[4];
		std::istringstream iss(check2);
		iss >> comp[1]; 
		iss >> comp[2];
		iss >> comp[3];
		if((xyzpdb->x(ii) - comp[1]) > 0.01){
			cerr <<" Warning check x pos atom "<<ii<<" got "<<comp[1]<<" and "<<xyzpdb->x(ii)<<" from string "<<LINE1<<endl;
		}
		if((xyzpdb->y(ii) - comp[2]) > 0.01){
			cerr <<" Warning check y pos atom "<<ii<<" got "<<comp[2]<<" and "<<xyzpdb->y(ii)<<" from string "<<LINE1<<endl;
		}

		if((xyzpdb->z(ii) - comp[3]) > 0.01){
			cerr <<" Warning check z pos atom "<<ii<<" got "<<comp[3]<<" and "<<xyzpdb->z(ii)<<" from string "<<LINE1<<endl;
		}
		ii++; //atom index advances from successful loading
		//cout << "atom added, type " << atom[CNT1].type << " pol "<< atom[CNT1].bondPOL << endl;
		//cout <<","<<LINE1.substr(31,7)<<","<<LINE1.substr(39,7)<<","<<LINE1.substr(47,7)<<endl;
	}
	infile.close();
	if(xyzpdb->n() == MAX){
		cerr <<"max number of atoms reached."<< endl;
	}
return xyzpdb;
}

XYZ_VEC_CRD * loadcrd_xyzvec(int MAX, string filename){
	XYZ_VEC_CRD * xyz = new XYZ_VEC_CRD;
	ifstream infile;
    infile.open (filename.c_str());
	if(infile.fail()==true){ 
		cerr << "Input File "<<filename<<" Not Found" << endl; 
		infile.close();
		throw 99;
		return xyz; 
	}
	//Log(1) <<" input file .xyz selected \n";  
	//read comment lines in xyz
	string LINE1;
	getline(infile,LINE1);
	int num_af = (int)StoF(LINE1); //read number of atoms in file from line1
	getline(infile,LINE1);
	xyz->set_comment(LINE1);
	for(int ii=1;ii<=MAX;ii++){
		getline(infile,LINE1);
		string FRAG;
		stringstream ss(LINE1);
		vector<string> tokens;
		int cc=1;
		double xx,yy,zz;
		string a2; 
		while( ss >> FRAG){
			tokens.push_back(FRAG);
			if(cc==1) a2=FRAG ; //atomtype is 1st 
			else if(cc==2) xx=StoF(FRAG) ; //x coord 
			else if(cc==3) yy=StoF(FRAG); //y coord
			else if(cc==4) zz=StoF(FRAG); //z coord
			cc++;
		}
		
		if(cc<5){  //not full line
			infile.close();
			if(xyz->n() != num_af){
				cerr <<"Number of atoms found "<<xyz->n()<<" is different from number specified in file "<<num_af<<endl;
			}
			return xyz; //error in file or line
		}
		xyz->add_pos(a2,xx,yy,zz);
		if(infile.eof()) return xyz; 
	}
infile.close();
	if(xyz->n() == MAX){
		cerr <<"max number of atoms reached."<< endl;
	}
return xyz;
}

XYZ_VEC_CRD * loadcrd_xyzchr_vec(int MAX, string filename){
	string FN = __FUNCTION__;
	XYZ_VEC_CRD * xyzc = new XYZ_VEC_CRD;
	ifstream infile;
    infile.open (filename.c_str());
	if(infile.fail()==true){ 
		cerr << "Input File "<<filename<<" Not Found" << endl; 
		infile.close();
		throw 99;
		return xyzc; 
	}
	Log(1) <<" input file .xyzchr selected \n";  
	//read comment lines in xyzchr
	string LINE1;
	getline(infile,LINE1);
	int num_af = (int)StoF(LINE1); //read number of atoms in file from line1
	if(num_af < MAX) MAX = num_af+1; //read 1 extra line to make sure crd lines stop
	getline(infile,LINE1);
	xyzc->set_comment(LINE1);
	int nzchr=0,nzpol=0;
	for(int ii=1;ii<=MAX;ii++){ //
		getline(infile,LINE1);
		string FRAG;
		stringstream ss(LINE1);
		vector<string> tokens;
		int cc=1;
		double xx,yy,zz,ch=0.0,pp=0.0;
		string a2,re=""; 
		while( ss >> FRAG){
			tokens.push_back(FRAG);
			if(cc==1) a2=FRAG ; //atomtype is 1st 
			else if(cc==2) xx=StoF(FRAG) ; //x coord 
			else if(cc==3) yy=StoF(FRAG); //y coord
			else if(cc==4) zz=StoF(FRAG); //z coord
			else if(cc==5) ch=StoF(FRAG); //z coord
			else if(cc==6) pp=StoF(FRAG); //z coord
			else if(cc==7) re=FRAG; //z coord
			cc++;
		}
		if(cc>1 && cc<5 ){
			Log()<<"Last line read ="<<LINE1<<"\n";
		}
		if(cc<5 ){  //not full crd line
			if(xyzc->n() != num_af){
				Log(-1,FN,1) <<"Last Line read ->"<<LINE1<<"\n";
			}
			Log(1)<<"Read stopped on line. "<<xyzc->n()<<" atoms loaded, "<<nzchr<<" with non-zero charge and "
				<<nzpol<<" with non-zero polarizability.\n";
			break;
		}
		if(cc<6) Log() << "Charge not found for atom # "<<ii<<" using 0.0\n";
		if(cc<6) Log(1) << "Polarizability not found for atom # "<<ii<<" using 0.0\n";
		if(cc<8) Log(2) << "Residue or MM label not found for atom #"<<ii<<"\n";
		if(ch != 0.0) nzchr++;
		if(pp != 0.0) nzpol++;
		xyzc->add_pcp(a2,xx,yy,zz,ch,pp,re);
		if(infile.eof()) break; 
	}
	if(xyzc->n() >  num_af){
		Log(-1,FN,1) <<"Number of atoms found is greater than number specified in file!\n";
	}
	if(xyzc->n() != num_af && xyzc->n() != MAX ){
		Log(-1,FN,1) <<"Number of atoms found "<<xyzc->n()<<" is different from number specified in file "<<num_af<<"\n";
		xyzc->cclear();
		infile.close();
		return xyzc;
	}
	if(xyzc->n() == MAX){
		cerr <<"max number of atoms "<< MAX <<" was reached."<< endl;
	}
	Log(1)<<xyzc->n()<<" atoms loaded, "<<nzchr<<" with non-zero charge and "<<nzpol
		<<" with non-zero polarizability.\n";
	//read general data if present
	if(!infile.eof()) getline(infile,LINE1);
	if(LINE1 == "--DATA--"){
		if(!infile.eof()){
			Log(1)<<"Additional coordinate data found: ";
			//read periodic conditions
			for(int ii=1; ii<=9; ii++){
				if(infile.eof()) break;
				getline(infile,LINE1);
				if(LINE1.length()<1) break;
				xyzc->setdat(ii,StoF(LINE1) );
				Log(1)<<xyzc->dat(ii)<<", ";
			}
			Log(1)<<"\n";
		}
	}
	else{ // may have to read through many skipped atoms
		while(!infile.eof()){
		    getline(infile,LINE1);
			if(LINE1 == "--DATA--"){
				if(!infile.eof()){
					Log(1)<<"Additional coordinate data found: ";
					//read periodic conditions
					for(int ii=1; ii<=9; ii++){
						if(infile.eof()) break;
						getline(infile,LINE1);
						if(LINE1.length()<1) break;
						xyzc->setdat(ii,StoF(LINE1) );
						Log(1)<<xyzc->dat(ii)<<", ";
					}
					Log(1)<<"\n";
				}
			}
		}
	}
	infile.close();
return xyzc;
}

XYZ_VEC_CRD * loadcrd_molvec(int MAX, string filename){
	string FN = __FUNCTION__;
	XYZ_VEC_CRD * xyz = new XYZ_VEC_CRD;
	ifstream infile;
    infile.open (filename.c_str());
	if(infile.fail()==true){ 
		Log(-1,FN,3) << "Input File "<<filename<<" Not Found.\n"; 
		infile.close();
		throw 99;
		return xyz; 
	}
	Log(1) <<" input file .xyz selected \n";  
	//read comment lines in xyz
	string LINE1;
	getline(infile,LINE1); getline(infile,LINE1); getline(infile,LINE1); getline(infile,LINE1);
	xyz->set_comment(LINE1);
	double atnum = (int)StoF( LINE1.substr(0,3) );
	Log(2) << "Read in "<<atnum<<" atomsz.\n";
	for(int ii=1;ii<=atnum;ii++){
		getline(infile,LINE1);
		string tt;
		double xx,yy,zz;
		xx = StoF(LINE1.substr(0,11));  
		yy = StoF(LINE1.substr(12,9));
		zz = StoF(LINE1.substr(22,9));
		tt = LINE1.substr(31,2);
		if(LINE1.length() < 40){  //not full line
			infile.close();
			if(xyz->n() < atnum){
				Log(-1,FN,1) <<"Number of atoms found "<<xyz->n()<<" is different from number specified in file "<<atnum<<"\n";
			}
			return xyz; //error in file or line
		}
		xyz->add_pos(tt,xx,yy,zz);
		if(infile.eof()){ 
			infile.close();
			if(xyz->n() < atnum){
				Log(-1,FN,1) <<"Number of atoms found "<<xyz->n()<<" is different from number specified in file "<<atnum<<"\n";
			}
			return xyz; //error in file or line		 
		}
	}
infile.close();
	if(xyz->n() == MAX){
		cerr <<"max number of atoms reached."<< endl;
	}
return xyz;
}

XYZ_VEC_CRD * loadcrd_amberRST(int Max, string topFN, string rstFN){
	string FN = __FUNCTION__;
	XYZ_VEC_CRD * xyz = new XYZ_VEC_CRD;
	//open top and rst files
	ifstream topf;
    topf.open (topFN.c_str());
	if(topf.fail()==true){ 
		Log(-1,FN,3) << "Input File "<<topFN<<" Not Found.\n"; 
		topf.close();
		throw 99;
		return xyz; 
	}
	ifstream rstf;
    rstf.open (rstFN.c_str());
	if(rstf.fail()==true){ 
		Log(-1,FN,3) << "Input File "<<rstFN<<" Not Found.\n"; 
		rstf.close();
		throw 99;
		return xyz; 
	}
	Log(1) <<" Amber format input file selected \n";  
	
	string LINE;
	int natom=0,nfind=0;
	
	while(!topf.eof()){
		getline(topf,LINE);
		if(LINE.find("%FLAG POINTERS")<2){ //get number of atoms to load
			getline(topf,LINE); //read format line
			getline(topf,LINE);
			natom = (int)StoF(LINE.substr(0,8) );
			break;
		}
	}
	Log(2)<<"Load data for "<<natom<<" atoms.\n";
	if(natom==0){
		Log(-1,FN,1)<<"Could not find number of atoms in PRMTOP file "<<topFN<<"\n";
		xyz->setn(0);
		return xyz;
	}
	while(!topf.eof()){
		getline(topf,LINE);
		if(LINE.find("%FLAG CHARGE")<2){ //get charges
			getline(topf,LINE); //read format line
			nfind=0;
			Log(2)<<"Reading Charges ";
			while(!topf.eof()){
				getline(topf,LINE);
				stringstream ssl(LINE);
				while(!ssl.eof()){
					nfind++;
					double chr;
					ssl >> chr;
					chr = chr/18.2223; //convert to atomic charge
					xyz->add_pcp("",0,0,0,chr,0,"");
					if(nfind<20) Log(2)<<xyz->c(nfind)<<",";
					if(nfind == natom) break;
				}
				if(nfind == natom) break;
			}
			break;
		}
	}
	if(natom != nfind){
		Log(-1,FN,1)<<"\nCould not find correct number of atom charges in PRMTOP file "<<topFN<<"\n";
		xyz->setn(0);
		return xyz;
	}	
	Log(1)<<"\nfound "<<nfind<<" atomic charges\n";
	nfind=0;
	while(!topf.eof()){
		getline(topf,LINE);
		if(LINE.find("%FLAG ATOMIC_NUMBER")<2){ //get charges
			getline(topf,LINE); //read format line
			while(!topf.eof()){
				getline(topf,LINE);
				stringstream ssl(LINE);
				while(!ssl.eof()){
					nfind++;
					int atomicn;
					ssl >> atomicn;
					string atype = AssignType(atomicn);
					xyz->sett(nfind,atype);
					if(nfind == natom) break;
				}
				if(nfind == natom) break;
			}
			break;
		}
	}
	Log(1)<<"\nfound "<<nfind<<" atomic numbers and converted them to atomic symbols\n";
	if(natom != nfind){
		Log(-1,FN,1)<<"Could not find correct number of atomic numbers in PRMTOP file "<<topFN<<"\n";
		xyz->setn(0);
		return xyz;
	}	
	topf.close();
	//read crd from rst file.
	getline(rstf,LINE); //read comment line
	streampos stpos = rstf.tellg();
	int rst_anum ; //how many sets of crd data to read
	getline(rstf,LINE);
	if(LINE.length() < 23){ //old file format
		rst_anum= (int)StoF( LINE.substr(0,5));
		if(rst_anum != natom){
			Log(-1,FN,1)<<"Different number of atoms in PRMTOP "<<natom<<" and RST file "<<rst_anum<<"?\n";
			xyz->setn(0);
			return xyz;	
		}
	}
	else{ //new file format
		rst_anum = natom;
		Log(0)<<"WARNING new file format for rst files. Cannot determine mismatch in rst and prmtop files!\n";
		//reset file
		rstf.seekg(stpos);
	}
	int idx=1,adx=1,lidx=0;
	while(!rstf.eof()){
		getline(rstf,LINE); 
		lidx++; //count number of lines	
		stringstream ssl(LINE);
		while(!ssl.eof()){	
			double point;
			ssl >> point;
			if(idx==1) xyz->setx(adx,point);  
			if(idx==2) xyz->sety(adx,point);
			if(idx==3) xyz->setz(adx,point);
			idx++;
			if(idx > 3){ idx=1; adx++; }
			if(adx > natom) break;
		}
		if(adx > natom) break;
	}
	stpos = rstf.tellg(); //get current position
	getline(rstf,LINE); //read the next line
	getline(rstf,LINE);  //read the next line which might bring file to end
	if(natom < 5){ Log(-1,FN,1)<<"Less than 5 atoms in rst file?? cannot properly read box parameters.\n"; }
	int readl = 2;
	if(!rstf.eof()){    //velocity info found
		Log(1)<<"Velocities found. Skipping "<<lidx<<" lines.\n";
		while(readl < lidx){
			getline(rstf,LINE);
			readl++;
		}
	}
	else{ //no velocities
		Log(1)<<"No Velocities found.\n";
		rstf.clear();
		rstf.seekg(stpos); //reset file position to read periodic boundary data
	}
	getline(rstf,LINE); 
	Log(2)<<"Last line:"<<LINE<<"\n";
	stringstream ssl2(LINE);
	Log(1)<<"Reading in Periodic boundary data...\n";
	int idx2=1;
	while(!ssl2.eof()){	
		double point;
		ssl2 >> point;
		Log(1)<<idx2<<"="<<point<<"\n";
		xyz->setdat(idx2,point);
		idx2++;
		
		if(idx2 > 9){
			Log(-1,FN,1)<<"rst file with greater than 9 periodic conditions??";
			break;
		}
	}
	rstf.close();
return xyz;
}




int printcrd_xyzchr_vec(XYZ_VEC_CRD * xyzp_crd, string filename1, string comment1){
	ofstream outfile2 ( filename1.c_str() );
	int CNT6 = xyzp_crd->n();
	outfile2 << CNT6 << endl;
	outfile2 << comment1 << endl;
	outfile2 << fixed <<showpoint;
	//outfile2 << setprecision(5);
	char XP[14],YP[14],ZP[14],CH[14],PO[14];
	for( int xx=1; xx<= CNT6; xx++){ // print xyz coor
			sprintf(XP,"%12.6f",xyzp_crd->x(xx) );
			sprintf(YP,"%12.6f",xyzp_crd->y(xx) );
			sprintf(ZP,"%12.6f",xyzp_crd->z(xx) );
			sprintf(CH,"%12.6f",xyzp_crd->c(xx) );
			sprintf(PO,"%12.6f",xyzp_crd->p(xx) );
			//outfile2  << xyzp_crd->type(xx) << "\t"  << xyzp_crd->x(xx) << "\t" << xyzp_crd->y(xx) << "\t" << xyzp_crd->z(xx) << endl;
			outfile2 << xyzp_crd->type(xx) << "\t"  << XP << "\t"  << YP << "\t"  << ZP << "\t" \
				<< CH << "\t"  << PO << "\t"  << xyzp_crd->gres(xx) <<endl;
	}
	//write data
	outfile2 <<endl<<"--DATA--"<<endl; //data marker
	for(int xx=1;xx<=9;xx++){
		outfile2 << xyzp_crd->dat(xx)<<endl;
	}
	outfile2.close();
	Log() <<"Wrote "<<CNT6<<" to "<<filename1<<"\n";
return 0;
}

int printcrd_xyzvec(XYZ_VEC_CRD * xyzp_crd, string filename1, string comment1){
	ofstream outfile2 ( filename1.c_str() );
	int CNT6 = xyzp_crd->n();
	outfile2 << CNT6 << endl;
	outfile2 << comment1 << endl;
	outfile2 << fixed <<showpoint;
	//outfile2 << setprecision(5);
	char XP[14],YP[14],ZP[14];
	for( int xx=1; xx<= CNT6; xx++){ // print xyz coor
			sprintf(XP,"%12.6f",xyzp_crd->x(xx) );
			sprintf(YP,"%12.6f",xyzp_crd->y(xx) );
			sprintf(ZP,"%12.6f",xyzp_crd->z(xx) );
			//outfile2  << xyzp_crd->type(xx) << "\t"  << xyzp_crd->x(xx) << "\t" << xyzp_crd->y(xx) << "\t" << xyzp_crd->z(xx) << endl;
			outfile2 << xyzp_crd->type(xx) << "\t"  << XP << "\t"  << YP << "\t"  << ZP <<endl;
	}
	outfile2.close();
	Log() <<"Wrote "<<CNT6<<" to "<<filename1<<"\n";
return 0;
}

int printcrd_xyzvec_precise(XYZ_VEC_CRD * xyzp_crd, string filename1, string comment1){
	ofstream outfile2 ( filename1.c_str() );
	int CNT6 = xyzp_crd->n();
	outfile2 << CNT6 << endl;
	outfile2 << comment1 << endl;
	outfile2 << fixed <<showpoint;
	//outfile2 << setprecision(5);
	char XP[25],YP[25],ZP[25];
	for( int xx=1; xx<= CNT6; xx++){ // print xyz coor
			sprintf(XP,"%24.18f",xyzp_crd->x(xx) );
			sprintf(YP,"%24.18f",xyzp_crd->y(xx) );
			sprintf(ZP,"%24.18f",xyzp_crd->z(xx) );
			//outfile2  << xyzp_crd->type(xx) << "\t"  << xyzp_crd->x(xx) << "\t" << xyzp_crd->y(xx) << "\t" << xyzp_crd->z(xx) << endl;
			outfile2 << xyzp_crd->type(xx) << "\t"  << XP << "\t"  << YP << "\t"  << ZP <<endl;
	}
	outfile2.close();
	Log() <<"Wrote "<<CNT6<<" to "<<filename1<<"\n";
return 0;
}

int printcrd_xyzvec_range(XYZ_VEC_CRD * xyzp_crd,int str,int fin, string filename1, string comment1){
	//prints out to file custom range
	ofstream outfile2 ( filename1.c_str() );
	int CNT6 = xyzp_crd->n();
	if(fin > CNT6){ 
		Log(-1) <<"Index "<<fin<<" is greater than number of atoms "<<CNT6<<"\n";
		return 0;
	}
	if(str < 1){ 
		Log(-1) <<"Start Index="<<str<<" is less than zero \n";
		return 0;
	}
	if(fin < str){ 
		Log(-1) <<"Index "<<str<<" is greater than end index "<<fin<<"\n";
		return 0;
	}	
	outfile2 << fin-str+1 << endl;
	outfile2 << comment1 << endl;
	outfile2 << fixed <<showpoint;
	//outfile2 << setprecision(5);
	char XP[14],YP[14],ZP[14];
	for( int xx=str; xx<= fin; xx++){ // print xyz coor
			sprintf(XP,"%12.6f",xyzp_crd->x(xx) );
			sprintf(YP,"%12.6f",xyzp_crd->y(xx) );
			sprintf(ZP,"%12.6f",xyzp_crd->z(xx) );
			//outfile2  << xyzp_crd->type(xx) << "\t"  << xyzp_crd->x(xx) << "\t" << xyzp_crd->y(xx) << "\t" << xyzp_crd->z(xx) << endl;
			outfile2 << xyzp_crd->type(xx) << "\t"  << XP << "\t"  << YP << "\t"  << ZP <<endl;
	}
	outfile2.close();
	Log() <<"Wrote "<<fin-str+1<<" to "<<filename1<<"\n";
return 0;
}

int type2atomN(string type){
     if ( type == "H" ) return 1;
else if ( type == "He" ) return 2;
else if ( type == "Li" ) return 3;
else if ( type == "Be" ) return 4;
else if ( type == "B" ) return 5;
else if ( type == "C" ) return 6;
else if ( type == "N" ) return 7;
else if ( type == "O" ) return 8;
else if ( type == "F" ) return 9;
else if ( type == "Ne" ) return 10;
else if ( type == "Na" ) return 11;
else if ( type == "Mg" ) return 12;
else if ( type == "Al" ) return 13;
else if ( type == "Si" ) return 14;
else if ( type == "P" ) return 15;
else if ( type == "S" ) return 16;
else if ( type == "Cl" ) return 17;
else if ( type == "Ar" ) return 18;
else if ( type == "K" ) return 19;
else if ( type == "Ca" ) return 20;
else if ( type == "Sc" ) return 21;
else if ( type == "Ti" ) return 22;
else if ( type == "V" ) return 23;
else if ( type == "Cr" ) return 24;
else if ( type == "Mn" ) return 25;
else if ( type == "Fe" ) return 26;
else if ( type == "Co" ) return 27;
else if ( type == "Ni" ) return 28;
else if ( type == "Cu" ) return 29;
else if ( type == "Zn" ) return 30;
else if ( type == "Ga" ) return 31;
else if ( type == "Ge" ) return 32;
else if ( type == "As" ) return 33;
else if ( type == "Se" ) return 34;
else if ( type == "Br" ) return 35;
else if ( type == "Kr" ) return 36;
else if ( type == "Rb" ) return 37;
else if ( type == "Sr" ) return 38;
else if ( type == "Y" ) return 39;
else if ( type == "Zr" ) return 40;
else if ( type == "Nb" ) return 41;
else if ( type == "Mo" ) return 42;
else if ( type == "Tc" ) return 43;
else if ( type == "Ru" ) return 44;
else if ( type == "Rh" ) return 45;
else if ( type == "Pd" ) return 46;
else if ( type == "Ag" ) return 47;
else if ( type == "Cd" ) return 48;
else if ( type == "In" ) return 49;
else if ( type == "Sn" ) return 50;
else if ( type == "Sb" ) return 51;
else if ( type == "Te" ) return 52;
else if ( type == "I" ) return 53;
else if ( type == "Xe" ) return 54;
else if ( type == "Cs" ) return 55;
else if ( type == "Ba" ) return 56;
else if ( type == "La" ) return 57;
else if ( type == "Ce" ) return 58;
else if ( type == "Pr" ) return 59;
else if ( type == "Nd" ) return 60;
else if ( type == "Pm" ) return 61;
else if ( type == "Sm" ) return 62;
else if ( type == "Eu" ) return 63;
else if ( type == "Gd" ) return 64;
else if ( type == "Tb" ) return 65;
else if ( type == "Dy" ) return 66;
else if ( type == "Ho" ) return 67;
else if ( type == "Er" ) return 68;
else if ( type == "Tm" ) return 69;
else if ( type == "Yb" ) return 70;
else if ( type == "Lu" ) return 71;
else if ( type == "Hf" ) return 72;
else if ( type == "Ta" ) return 73;
else if ( type == "W" ) return 74;
else if ( type == "Re" ) return 75;
else if ( type == "Os" ) return 76;
else if ( type == "Ir" ) return 77;
else if ( type == "Pt" ) return 78;
else if ( type == "Au" ) return 79;
else if ( type == "Hg" ) return 80;
else if ( type == "Tl" ) return 81;
else if ( type == "Pb" ) return 82;
else if ( type == "Bi" ) return 83;
else if ( type == "Po" ) return 84;
else if ( type == "At" ) return 85;
else if ( type == "Rn" ) return 86;
else if ( type == "Fr" ) return 87;
else if ( type == "Ra" ) return 88;
else if ( type == "Ac" ) return 89;
else if ( type == "Th" ) return 90;
else if ( type == "Pa" ) return 91;
else if ( type == "U" ) return 92;
else if ( type == "Np" ) return 93;
else if ( type == "Pu" ) return 94;
else if ( type == "Am" ) return 95;
else if ( type == "Cm" ) return 96;
else if ( type == "Bk" ) return 97;
else if ( type == "Cf" ) return 98;
else if ( type == "Es" ) return 99;
else if ( type == "Fm" ) return 100;
else if ( type == "Md" ) return 101;
else if ( type == "No" ) return 102;
else if ( type == "Lr" ) return 103;
else if ( type == "Rf" ) return 104;
else if ( type == "Db" ) return 105;
else if ( type == "Sg" ) return 106;
else if ( type == "Bh" ) return 107;
else if ( type == "Hs" ) return 108;
else if ( type == "Mt" ) return 109;
else if ( type == "Ds" ) return 110;
else if ( type == "Rg" ) return 111;
else if ( type == "Cn" ) return 112;
else if ( type == "Uut" ) return 113;
else if ( type == "Fl" ) return 114;
else if ( type == "Uup" ) return 115;
else if ( type == "Lv" ) return 116;
else if ( type == "Uus" ) return 117;
else if ( type == "Uuo" ) return 118;
    else cerr << "WARNING!!!!   Unidentified Atom Type: "<<type<< " Consider adding atomtypes" <<endl;
return 0;
}

string AssignType(int AtomN){
     if(AtomN == 1   ) return "H " ; //   Hydrogen
else if(AtomN == 2   ) return "He" ; //    Helium
else if(AtomN == 3   ) return "Li" ; //    Lithium
else if(AtomN == 4   ) return "Be" ; //    Beryllium
else if(AtomN == 5   ) return "B " ; //   Boron
else if(AtomN == 6   ) return "C " ; //   Carbon
else if(AtomN == 7   ) return "N " ; //   Nitrogen
else if(AtomN == 8   ) return "O " ; //   Oxygen
else if(AtomN == 9   ) return "F " ; //   Fluorine
else if(AtomN == 10  ) return "Ne" ; //    Neon
else if(AtomN == 11  ) return "Na" ; //    Sodium
else if(AtomN == 12  ) return "Mg" ; //    Magnesium
else if(AtomN == 13  ) return "Al" ; //    Aluminum, Aluminium
else if(AtomN == 14  ) return "Si" ; //    Silicon
else if(AtomN == 15  ) return "P " ; //   Phosphorus
else if(AtomN == 16  ) return "S " ; //   Sulfur
else if(AtomN == 17  ) return "Cl" ; //    Chlorine
else if(AtomN == 18  ) return "Ar" ; //    Argon
else if(AtomN == 19  ) return "K " ; //   Potassium
else if(AtomN == 20  ) return "Ca" ; //    Calcium
else if(AtomN == 21  ) return "Sc" ; //    Scandium
else if(AtomN == 22  ) return "Ti" ; //    Titanium
else if(AtomN == 23  ) return "V " ; //   Vanadium
else if(AtomN == 24  ) return "Cr" ; //    Chromium
else if(AtomN == 25  ) return "Mn" ; //    Manganese
else if(AtomN == 26  ) return "Fe" ; //    Iron
else if(AtomN == 27  ) return "Co" ; //    Cobalt
else if(AtomN == 28  ) return "Ni" ; //    Nickel
else if(AtomN == 29  ) return "Cu" ; //    Copper
else if(AtomN == 30  ) return "Zn" ; //    Zinc
else if(AtomN == 31  ) return "Ga" ; //    Gallium
else if(AtomN == 32  ) return "Ge" ; //    Germanium
else if(AtomN == 33  ) return "As" ; //    Arsenic
else if(AtomN == 34  ) return "Se" ; //    Selenium
else if(AtomN == 35  ) return "Br" ; //    Bromine
else if(AtomN == 36  ) return "Kr" ; //    Krypton
else if(AtomN == 37  ) return "Rb" ; //    Rubidium
else if(AtomN == 38  ) return "Sr" ; //    Strontium
else if(AtomN == 39  ) return "Y " ; //   Yttrium
else if(AtomN == 40  ) return "Zr" ; //    Zirconium
else if(AtomN == 41  ) return "Nb" ; //    Niobium
else if(AtomN == 42  ) return "Mo" ; //    Molybdenum
else if(AtomN == 43  ) return "Tc" ; //    Technetium
else if(AtomN == 44  ) return "Ru" ; //    Ruthenium
else if(AtomN == 45  ) return "Rh" ; //    Rhodium
else if(AtomN == 46  ) return "Pd" ; //    Palladium
else if(AtomN == 47  ) return "Ag" ; //    Silver
else if(AtomN == 48  ) return "Cd" ; //    Cadmium
else if(AtomN == 49  ) return "In" ; //    Indium
else if(AtomN == 50  ) return "Sn" ; //    Tin
else if(AtomN == 51  ) return "Sb" ; //    Antimony
else if(AtomN == 52  ) return "Te" ; //    Tellurium
else if(AtomN == 53  ) return "I " ; //   Iodine
else if(AtomN == 54  ) return "Xe" ; //    Xenon
else if(AtomN == 55  ) return "Cs" ; //    Cesium
else if(AtomN == 56  ) return "Ba" ; //    Barium
else if(AtomN == 57  ) return "La" ; //    Lanthanum
else if(AtomN == 58  ) return "Ce" ; //    Cerium
else if(AtomN == 59  ) return "Pr" ; //    Praseodymium
else if(AtomN == 60  ) return "Nd" ; //    Neodymium
else if(AtomN == 61  ) return "Pm" ; //    Promethium
else if(AtomN == 62  ) return "Sm" ; //    Samarium
else if(AtomN == 63  ) return "Eu" ; //    Europium
else if(AtomN == 64  ) return "Gd" ; //    Gadolinium
else if(AtomN == 65  ) return "Tb" ; //    Terbium
else if(AtomN == 66  ) return "Dy" ; //    Dysprosium
else if(AtomN == 67  ) return "Ho" ; //    Holmium
else if(AtomN == 68  ) return "Er" ; //    Erbium
else if(AtomN == 69  ) return "Tm" ; //    Thulium
else if(AtomN == 70  ) return "Yb" ; //    Ytterbium
else if(AtomN == 71  ) return "Lu" ; //    Lutetium
else if(AtomN == 72  ) return "Hf" ; //    Hafnium
else if(AtomN == 73  ) return "Ta" ; //    Tantalum
else if(AtomN == 74  ) return "W " ; //   Tungsten
else if(AtomN == 75  ) return "Re" ; //    Rhenium
else if(AtomN == 76  ) return "Os" ; //    Osmium
else if(AtomN == 77  ) return "Ir" ; //    Iridium
else if(AtomN == 78  ) return "Pt" ; //    Platinum
else if(AtomN == 79  ) return "Au" ; //    Gold
else if(AtomN == 80  ) return "Hg" ; //    Mercury
else if(AtomN == 81  ) return "Tl" ; //    Thallium
else if(AtomN == 82  ) return "Pb" ; //    Lead
else if(AtomN == 83  ) return "Bi" ; //    Bismuth
else if(AtomN == 84  ) return "Po" ; //    Polonium
else if(AtomN == 85  ) return "At" ; //    Astatine
else if(AtomN == 86  ) return "Rn" ; //    Radon
else if(AtomN == 87  ) return "Fr" ; //    Francium
else if(AtomN == 88  ) return "Ra" ; //    Radium
else if(AtomN == 89  ) return "Ac" ; //    Actinium
else if(AtomN == 90  ) return "Th" ; //    Thorium
else if(AtomN == 91  ) return "Pa" ; //    Protactinium
else if(AtomN == 92  ) return "U " ; //   Uranium
else if(AtomN == 93  ) return "Np" ; //    Neptunium
else if(AtomN == 94  ) return "Pu" ; //    Plutonium
else if(AtomN == 95  ) return "Am" ; //    Americium
else if(AtomN == 96  ) return "Cm" ; //    Curium
else if(AtomN == 97  ) return "Bk" ; //    Berkelium
else if(AtomN == 98  ) return "Cf" ; //    Californium
else if(AtomN == 99  ) return "Es" ; //    Einsteinium
else if(AtomN == 100 ) return "Fm "; //   Fermium
else if(AtomN == 101 ) return "Md "; //   Mendelevium
else if(AtomN == 102 ) return "No "; //   Nobelium
else if(AtomN == 103 ) return "Lr "; //   Lawrencium
else if(AtomN == 104 ) return "Rf "; //   Rutherfordium
else if(AtomN == 105 ) return "Db "; //   Dubnium
else if(AtomN == 106 ) return "Sg "; //   Seaborgium
else if(AtomN == 107 ) return "Bh "; //   Bohrium
else if(AtomN == 108 ) return "Hs "; //   Hassium
else if(AtomN == 109 ) return "Mt "; //   Meitnerium
else if(AtomN == 110 ) return "Ds "; //   Darmstadtium
else if(AtomN == 111 ) return "Rg "; //   Roentgenium
else if(AtomN == 112 ) return "Cn "; //   Copernicium
else if(AtomN == 113 ) return "Uut"; //   Ununtrium
else if(AtomN == 114 ) return "Fl "; //   Flerovium
else if(AtomN == 115 ) return "Uup"; //   Ununpentium
else if(AtomN == 116 ) return "Lv "; //  Livermorium
else if(AtomN == 117 ) return "Uus"; //   Ununseptium
else if(AtomN == 118 ) return "Uuo"; //   Ununoctium
    else cout << "WARNING!!!!   Unidentified Atom Type, Atomic Number "<<AtomN<< " Consider adding atomtypes" <<endl;
return "Du";
}

double AssignMass(string type){ // return Isotopic mass from atom type
if ( type == "H" ) return 1.0079;
else if ( type == "He" ) return 4.0026;
else if ( type == "Li" ) return 6.941;
else if ( type == "Be" ) return 9.0122;
else if ( type == "B" ) return 10.811;
else if ( type == "C" ) return 12.0107;
else if ( type == "N" ) return 14.0067;
else if ( type == "O" ) return 15.9994;
else if ( type == "F" ) return 18.9984;
else if ( type == "Ne" ) return 20.1797;
else if ( type == "Na" ) return 22.9897;
else if ( type == "Mg" ) return 24.305;
else if ( type == "Al" ) return 26.9815;
else if ( type == "Si" ) return 28.0855;
else if ( type == "P" ) return 30.9738;
else if ( type == "S" ) return 32.065;
else if ( type == "Cl" ) return 35.453;
else if ( type == "Ar" ) return 39.948;
else if ( type == "K" ) return 39.0983;
else if ( type == "Ca" ) return 40.078;
else if ( type == "Sc" ) return 44.9559;
else if ( type == "Ti" ) return 47.867;
else if ( type == "V" ) return 50.9415;
else if ( type == "Cr" ) return 51.9961;
else if ( type == "Mn" ) return 54.938;
else if ( type == "Fe" ) return 55.845;
else if ( type == "Co" ) return 58.9332;
else if ( type == "Ni" ) return 58.6934;
else if ( type == "Cu" ) return 63.546;
else if ( type == "Zn" ) return 65.39;
else if ( type == "Ga" ) return 69.723;
else if ( type == "Ge" ) return 72.64;
else if ( type == "As" ) return 74.9216;
else if ( type == "Se" ) return 78.96;
else if ( type == "Br" ) return 79.904;
else if ( type == "Kr" ) return 83.8;
else if ( type == "Rb" ) return 85.4678;
else if ( type == "Sr" ) return 87.62;
else if ( type == "Y" ) return 88.9059;
else if ( type == "Zr" ) return 91.224;
else if ( type == "Nb" ) return 92.9064;
else if ( type == "Mo" ) return 95.94;
else if ( type == "Tc" ) return 98;
else if ( type == "Ru" ) return 101.07;
else if ( type == "Rh" ) return 102.9055;
else if ( type == "Pd" ) return 106.42;
else if ( type == "Ag" ) return 107.8682;
else if ( type == "Cd" ) return 112.411;
else if ( type == "In" ) return 114.818;
else if ( type == "Sn" ) return 118.71;
else if ( type == "Sb" ) return 121.76;
else if ( type == "Te" ) return 127.6;
else if ( type == "I" ) return 126.9045;
else if ( type == "Xe" ) return 131.293;
else if ( type == "Cs" ) return 132.9055;
else if ( type == "Ba" ) return 137.327;
else if ( type == "La" ) return 138.9055;
else if ( type == "Ce" ) return 140.116;
else if ( type == "Pr" ) return 140.9077;
else if ( type == "Nd" ) return 144.24;
else if ( type == "Pm" ) return 145;
else if ( type == "Sm" ) return 150.36;
else if ( type == "Eu" ) return 151.964;
else if ( type == "Gd" ) return 157.25;
else if ( type == "Tb" ) return 158.9253;
else if ( type == "Dy" ) return 162.5;
else if ( type == "Ho" ) return 164.9303;
else if ( type == "Er" ) return 167.259;
else if ( type == "Tm" ) return 168.9342;
else if ( type == "Yb" ) return 173.04;
else if ( type == "Lu" ) return 174.967;
else if ( type == "Hf" ) return 178.49;
else if ( type == "Ta" ) return 180.9479;
else if ( type == "W" ) return 183.84;
else if ( type == "Re" ) return 186.207;
else if ( type == "Os" ) return 190.23;
else if ( type == "Ir" ) return 192.217;
else if ( type == "Pt" ) return 195.078;
else if ( type == "Au" ) return 196.9665;
else if ( type == "Hg" ) return 200.59;
else if ( type == "Tl" ) return 204.3833;
else if ( type == "Pb" ) return 207.2;
else if ( type == "Bi" ) return 208.9804;
else if ( type == "Po" ) return 209;
else if ( type == "At" ) return 210;
else if ( type == "Rn" ) return 222;
else if ( type == "Fr" ) return 223;
else if ( type == "Ra" ) return 226;
else if ( type == "Ac" ) return 227;
else if ( type == "Th" ) return 232.0381;
else if ( type == "Pa" ) return 231.0359;
else if ( type == "U" ) return 238.0289;
else if ( type == "Np" ) return 237;
else if ( type == "Pu" ) return 244;
else if ( type == "Am" ) return 243;
else if ( type == "Cm" ) return 247;
else if ( type == "Bk" ) return 247;
else if ( type == "Cf" ) return 251;
else if ( type == "Es" ) return 252;
else if ( type == "Fm" ) return 257;
else if ( type == "Md" ) return 258;
else if ( type == "No" ) return 259;
else if ( type == "Lr" ) return 262;
else if ( type == "Rf" ) return 261;
else if ( type == "Db" ) return 262;
else if ( type == "Sg" ) return 266;
else if ( type == "Bh" ) return 264;
else if ( type == "Hs" ) return 277;
else if ( type == "Mt" ) return 268;
 else cerr << "WARNING!!!!   Unidentified Atom Type  "<<type<< " Not found. Consider adding atomtypes" <<endl;

return 0;
}

inline std::string trimstr(std::string& str){
str.erase(0, str.find_first_not_of(' '));       //prefixing spaces
str.erase(str.find_last_not_of(' ')+1);         //surfixing spaces
return str;
}
