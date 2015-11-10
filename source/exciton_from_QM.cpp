#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <math.h>
#include <vector>
using namespace std;
 
#include "class_log.h" //custom print message class
#include "load_crdf.h"  // load xyz_crd class
#include "xyz_manip.h" // manipulate  crd 
#include "xyz_analysis.h" // 
const double O16M=16.00; //mass of oxygen
const double C12M=12.00; //mass of Carbon
const double warn_disp=0.85; // warn if selected atoms summed displacements is below this value

string COMMENT2, HELPLINE="Program written to analyze the exciton coupling from a gaussian freq=vcd output file. \n\
Usage: < -i {g09.out} g09 Freq output file> <-n1 (normal mode 1)> <-n2 (normal mode 2)> <-h help>\n\
	   < -g1 (group 1, Cabon Atom) (group 1, Oxygen atom)> < -g2 (group 2, Cabon Atom) (group 2, Oxygen atom)>\n\
	   < -l (integer 0-3, loudness) > < -p just print moving atoms and quit> < -a (integer, atom number) print atom >\n\
Use of the '-a' option will print additional atoms displacements, and can be used for multiple atoms,\n\
	for example -a 3 -a 4 -a 7   will add atoms 3,4, and 7 to the normal mode print list. \n\
TIP: use -l 2 to get information in table format.\n\n\
Please cite \"Determination of the Absolute Configurations Using\n\
   Exciton Chirality Method for Vibrational Circular Dichroism:\n\
   Right Answers for the Wrong Reasons?, C. L. Covington, V. P. Nicu,\n\
   P.L.Polavarapu,  J. Phys. Chem. A, 2015, 119 (42), pp 10589-10601\"\n";
//	   < --detect use automatic detection of normal mode and groups>\n";
	   
string PN="exciton_from_QM";

string mod_date="6-16-2015", cmp_date=__DATE__;
int loudness = 1; //set with -l  0 very quite, 1 normal, 2 noisy, 3 debug , -1 error
#define FN __FUNCTION__
// Changed on 12-15-2014 to XYZ_VEC functions
 
int main (int argc, char *argv[])
{     
    if (argc < 2)
    {// If the user didn't provide a filename command line argument, print an error and exit.
        cerr << PN<<" requires arguments. Use -h for help\n";
        return 1;
    } 
    string INFILE="g09.out"; //default file in and out 
    int i=1,nm1=0,nm2=0,g1a1=0,g1a2=0,g2a1=0,g2a2=0;
	bool early_quit=false;
	vector<int> print_list;
    while (i + 1 <= argc){ // Check that we haven't finished parsing already
        string ARG=argv[i];
        if (ARG.find("-i")<2) {     INFILE=argv[i+1]; } 
        else if (ARG.find("-n1")<2) { nm1= (int)StoF(argv[i+1]); } 
		else if (ARG.find("-n2")<2) { nm2= (int)StoF(argv[i+1]); } 
		else if (ARG.find("-g1")<2) { g1a1= (int)StoF(argv[i+1]); g1a2= (int)StoF(argv[i+2]);  } 
		else if (ARG.find("-g2")<2) { g2a1= (int)StoF(argv[i+1]); g2a2= (int)StoF(argv[i+2]);  } 
		else if (ARG.find("-l")<1) { loudness= (int)StoF(argv[i+1]); }
		else if (ARG.find("-h")<1) { cout << HELPLINE << endl; return 0; }
		else if (ARG.find("-p")<1) { early_quit=true; }
		else if (ARG.find("-a")<2) { print_list.push_back( (int)StoF(argv[i+1]) ) ; } 
        i++;
    }
    cout <<PN<< " by Cody Covington. Last modified "<< mod_date<<". Use -h for help." << endl<<"File:"<<INFILE<<endl;
    string STR;
	XYZ_VEC_CRD crd;
	ifstream infile;
    infile.open (INFILE.c_str());
	 if(infile.fail()==true){ cout << "Input File "<<INFILE<<" Not Found" << endl; return 1; }
	// load crd from g09 file
	long int crdpos=0; //position of coordinate information. 
	while(!infile.eof()){
		getline(infile,STR);
		if(STR.find("     Standard orientation:")!=string::npos){
			crdpos = infile.tellg();
			Log(2)<<"found set of xyz crd at "<<crdpos<<" bytes.\n";
		}
		if(STR.find(" Frequencies --")!=string::npos){
			break;
		}
	}
	if(infile.eof()){
		if(crdpos > 0 ) infile.clear(); // clear eof flag
		else{
			cerr <<"Error in g09 output file. No crd information found?\n";
			return 1;
		}
	}
	infile.seekg(crdpos,infile.beg); // goto last crd section in file 
	Log(2)<<"Using crd at "<<crdpos<<" bytes.\n";
	// get xyz crd
			getline(infile,STR); getline(infile,STR); 
			getline(infile,STR); getline(infile,STR);//read 4 lines 
			while(!infile.eof()){
				getline(infile,STR);  //cout<<STRING<<endl; 
				if(STR.find("-----")!=string::npos){  goto FOUNDCRD; } 
				stringstream ss(STR); double comp[7];
				ss >> comp[1];
				ss >> comp[2];
				ss >> comp[3];
				ss >> comp[4];
				ss >> comp[5];
				ss >> comp[6];
				double xx,yy,zz;
				int an  = StoF(STR.substr(14,4));  //atomic number
				string type = AssignType(an);
				xx = StoF(STR.substr(35,48));  if( fabs(xx - comp[4]) > 0.005){ cerr<<"Check x coordinates?\n";}
				yy = StoF(STR.substr(49,60));  if( fabs(yy - comp[5]) > 0.005){ cerr<<"Check y coordinates?\n";}
				zz = StoF(STR.substr(61,71));  if( fabs(zz - comp[6]) > 0.005){ cerr<<"Check z coordinates?\n";}
				crd.add_pcp(type,xx,yy,zz,0,0,"");
				//cout <<"Found crd :"<<type<<" "<<xx<<" "<<yy<<" "<<zz<<"\n";
			}

	FOUNDCRD:
	int num=crd.n();
	cout<<num<<" cartesian coordinates found.\n";
	if(num<4){
		cerr<<"error in cartesian crd. only loaded "<<num<<" atoms!\n";
		infile.close();
		return 1;
	}
	// get freq and dipole str
	vector<double> Freq; // frequency
	vector<double> DStr; // Dipole Strength
	vector<double> RStr; // Rotational Strength
	int l=0;
	double nm1_gdisp=0; //displacement maginitudes of selected atoms in normal mode1
	double nm2_gdisp=0; //displacement maginitudes of selected atoms in normal mode2
	cout<<"\nCarbonyl Region Normal Modes (1600-1950 cm-1):\n";
	Torder1 test_bv1 = getAP(&crd,g1a2) - getAP(&crd,g1a1); //bond vector for group 1, point from carbon to oxygen
	Torder1 test_bv2 = getAP(&crd,g2a2) - getAP(&crd,g2a1); //bond vector for group 2, point from carbon to oxygen
	int nm1_symS=1, nm2_symS=1; //determines if given modes are symmetric or antisymmetric
	while(!infile.eof()){
		getline(infile,STR);
		if(STR.find(" Frequencies --")!=string::npos){
			// get frequencies
			Freq.push_back( StoF(STR.substr(16,14)));  //cout << STR.substr(16,14)<<"\n";
			Freq.push_back( StoF(STR.substr(39,14)));  //cout << STR.substr(39,14)<<"\n";
			Freq.push_back( StoF(STR.substr(62,14)));  //cout << STR.substr(62,14)<<"\n";
			getline(infile,STR);getline(infile,STR);getline(infile,STR);getline(infile,STR); //load 4 lines
			DStr.push_back( StoF(STR.substr(16,14)));
			DStr.push_back( StoF(STR.substr(39,14)));
			DStr.push_back( StoF(STR.substr(62,14)));
			getline(infile,STR);
			RStr.push_back( StoF(STR.substr(16,14)));
			RStr.push_back( StoF(STR.substr(39,14)));
			RStr.push_back( StoF(STR.substr(62,14)));
			getline(infile,STR); 
			if(STR.find("  Atom  AN      X      Y      Z")==string::npos) getline(infile,STR); // get extra line in some versions
			long int nmpos = infile.tellg();;
			for(i=1;i<=3;i++){
				
				l++; 
				if( Freq.at(l-1)>1600 && Freq.at(l-1)<1950 ){ //go back and print  big movers
					infile.seekg(nmpos,infile.beg);
					cout <<"Mode:"<<l<<" Freq="<<Freq.at(l-1)<<" DS="<<DStr.at(l-1)<<" RS="<<RStr.at(l-1)<<"\n";
					int px=3,py=4,pz=5;
					if(i==2){ px=6; py=7; pz=8; }
					if(i==3){ px=9; py=10; pz=11; }
					double dX,dY,dZ;
					Torder1 disp_g1a1, disp_g1a2, disp_g2a1, disp_g2a2; // displaced bond vectors
					for(int n=1; n<=num; n++){
						getline(infile,STR);
						stringstream ssn(STR);
						double compn[15];
						for(int k=1;k<=14;k++){ 
							ssn >> compn[k] ;
							if(ssn.eof()) break;
						}
						dX=compn[px];
						dY=compn[py];
						dZ=compn[pz];
						double disp=sqrt( (dX*dX)+(dY*dY)+(dZ*dZ) );
						string atype1= AssignType(compn[2]);
						Torder1 atom_disp(dX,dY,dZ);
						if(n==g1a1 || n==g1a2 || n==g2a1 || n==g2a2){ // if a selected atom
							if(l==nm1) nm1_gdisp += disp;
							if(l==nm2) nm2_gdisp += disp;
							cout <<" *Selected Atoms: Atom="<<n<<" type="<<atype1<<" NM_disp=["<<dX<<","<<dY<<","<<dZ<<"]\tMag="<<disp<<"\n";
						}
						else if(disp > 0.25){ 
							cout <<"    Moving Atoms: Atom="<<n<<" type="<<atype1<<" NM_disp=["<<dX<<","<<dY<<","<<dZ<<"]\tMag="<<disp<<"\n";
						}
						else{
							for(unsigned int pp=0; pp < print_list.size() ;pp++){
								if( n == print_list.at(pp) ){
									cout <<"    *Watch Atoms: Atom="<<n<<" type="<<atype1<<" NM_disp=["<<dX<<","<<dY<<","<<dZ<<"]\tMag="<<disp<<"\n";
								}
							}
						}
						if(n==g1a1) disp_g1a1 = getAP(&crd,g1a1) + atom_disp;
						if(n==g1a2) disp_g1a2 = getAP(&crd,g1a2) + atom_disp;
						if(n==g2a1) disp_g2a1 = getAP(&crd,g2a1) + atom_disp;
						if(n==g2a2) disp_g2a2 = getAP(&crd,g2a2) + atom_disp;
						
					}
					int test_symS=1; // multiply by -1 if compressing
					Torder1 test_bv1_disp = disp_g1a1 - disp_g1a2;
					if(test_bv1_disp.mag() > test_bv1.mag() ) cout <<" *"<<g1a1<<"-"<<g1a2<<" Stretching.\n";
					else{ 
						cout <<" *"<<g1a1<<"-"<<g1a2<<" Compressing.\n"; 					
						test_symS *= (-1);
					}
					
					Torder1 test_bv2_disp = disp_g2a1 - disp_g2a2;
					if(test_bv2_disp.mag() > test_bv2.mag() ) cout <<" *"<<g2a1<<"-"<<g2a2<<" Stretching.\n";
					else{ 
						cout <<" *"<<g2a1<<"-"<<g2a2<<" Compressing.\n";					
						test_symS *= (-1);
					}
					if(test_symS==-1) cout <<" **Antisymmetric Stretch detected.\n";
					else              cout <<" **Symmetric Stretch detected.\n";
					
					if(l==nm1) nm1_symS=test_symS;
					if(l==nm2) nm2_symS=test_symS;
					
				}
			}
			
			// read normal mode displacements for autodetect???
			// not yet implemented
			
		}
		if(STR.find(" Mulliken charges:")!=string::npos){
			getline(infile,STR); // skip 1 line
			for(int n=1; n<=num; n++){
				if(infile.eof()) cerr<<"Could not get all mulliken charges?\n";
				getline(infile,STR);
				stringstream ss(STR); double comp;
				ss >> comp; 
				ss >> STR;
				ss >> comp;
				Log(3)<<" read Mullikn charge at atom "<<n<<" type "<<STR<<" to be "<<comp<<"\n";
				crd.setc(n,comp);
			}
			
		}
	}
	infile.close();
	if(early_quit==true) return 0;
	cout <<l<<" Normal modes loaded.\nAnalyze NMs "<<nm1<<" and "<<nm2<<"\n";
	cout <<"NM "<<nm1<<" atoms "<<g1a1<<" and "<<g1a2<<"\n";
	cout <<"NM "<<nm2<<" atoms "<<g2a1<<" and "<<g2a2<<"\n";
	if(nm1 > l || nm2 > l){
		cerr<<"Normal mode number exceeds number loaded!\n";
		return 1;
	}
	if(nm1==nm2){
		cerr<<"valid, different normal modes must be entered!\n";
		return 1;
	}
	if(g1a1 > num || g1a2 > num || g2a1 > num || g2a2 > num ){
		cerr<<" atom index exceeds number loaded!\n";
		return 1;
	}
	if(g1a1 < 1 || g1a2 < 1 || g2a1 < 1 || g2a2 < 1 ){
		cerr<<" atom index not valid!\n";
		return 1;
	}
// ----  sanity check
	if(nm1_gdisp < warn_disp) cerr<<"warning NM="<<nm1<<" total displacements for selected atoms is only "<<nm1_gdisp<<"\n";
	if(nm2_gdisp < warn_disp) cerr<<"warning NM="<<nm2<<" total displacements for selected atoms is only "<<nm2_gdisp<<"\n";
	cout <<"Normal Mode "<<nm1<<" all selected atoms displacements="<<nm1_gdisp<<"\n";
	cout <<"Normal Mode "<<nm2<<" all selected atoms displacements="<<nm2_gdisp<<"\n";
 //------------ do exciton coupling analysis---------------
	//get average Dipole strength printed in (10**-40 esu**2 cm**2)
	double D = (DStr.at(nm1-1) + DStr.at(nm2-1))/2 ; // average Dipole Str
	double Desu = D*(1E-40);
	double edtm = sqrt (D);
	  Log(2)<<"D="<<D<<", Desu="<<Desu<<", edtm="<<edtm<<"\n";
	Torder1 bv1 = getAP(&crd,g1a2) - getAP(&crd,g1a1); //bond vector for group 1, point from carbon to oxygen
	  Log(3)<<"Bond_Vec1="<<bv1;
	bv1.norm(); //make unit vector
	  Log(3)<<", Bond_Vec1_norm="<<bv1;
	Torder1 g1cm = ((O16M % getAP(&crd,g1a2)) + (C12M % getAP(&crd,g1a1)));  //group 1 center of mass
	g1cm %= (1/(C12M+O16M));												 //group 1 center of mass
	  Log(3)<<", Group1_COM="<<g1cm<<"\n";
	Torder1 bv2 = getAP(&crd,g2a2) - getAP(&crd,g2a1); //bond vector for group 2, point from carbon to oxygen
	  Log(3)<<"Bond_Vec2="<<bv2;
	bv2.norm(); //make unit vector	
	  Log(3)<<", Bond_Vec2_norm="<<bv2;
	Torder1 g2cm = ((O16M % getAP(&crd,g2a2)) + (C12M % getAP(&crd,g2a1)) ); //group 2 center of mass
	g2cm %= (1/(C12M+O16M));												 //group 2 center of mass
	  Log(3)<<", Group2_COM="<<g2cm<<"\n";
	Torder1 vecR12 = g2cm - g1cm;    // vector between center of masses in Angstroms
	  
	Torder1 vecR12cm = vecR12 % 1e-8;  // vector between center of masses in cm
	double R12 = vecR12.mag(); // Distance in Angstroms
	double R12cm = R12/(1E8);  // in cm
	  Log(3)<<"R_12="<<vecR12<<", R_12_Magnitude="<<R12<<"\n";
	//
	double V12_t1 = dotp(bv1,bv2)*Desu/pow(R12cm,3);
	  Log(3)<<"Dot Prod bv1*bv2="<< dotp(bv1,bv2) <<"\n";
	double DOT_PROD = acos(dotp(bv1,bv2))*180/3.141592;
	cout <<"R12="<<R12<<" edtm="<<edtm<<"\n";
	double V12_t2 = dotp( bv1,vecR12cm) * dotp(bv2,vecR12cm);
	  Log(3)<<"mu_1*R12 = "<< dotp( bv1,vecR12cm) <<" mu_2*R12 = "<< dotp(bv2,vecR12cm)<<"\n";
	V12_t2 *= (3*Desu/pow(R12cm,5)); 
	cout <<"V12 = "<< V12_t1 <<" - "<<V12_t2<<"\n";
	double V12 = V12_t1 - V12_t2;
	cout <<"V12 = "<< V12 <<" erg \n"; // 1 erg = 1 esu**2 cm**2
	// convert to wavenumbers
	//   h = 6.6260755e-27 erg s
	//
	V12 *= 1/(29979245800*(6.6260755e-27));
	cout <<" V12 = "<<V12<<" cm-1\n";
	//  ------------- double check V12 using actually displaced bonds and mulliken charges.
		// Energy = const* Sum Q1*Q2/R12  -- ignore energy between atoms in same group.

	double Eexact[5];
	double facG1,facG2;
	for(int ql=0; ql<=4; ql++){
		Eexact[ql]=0;
		Torder1 pG1A1 = getAP(&crd,g1a1) ;
		Torder1 pG1A2 = getAP(&crd,g1a2) ;
		Torder1 pG2A1 = getAP(&crd,g2a1) ;
		Torder1 pG2A2 = getAP(&crd,g2a2) ;
		if(ql==0){ facG1=0.0; facG2=0.0;       Log(2)<<"E_EquilC ";   } // equilibrium crd
		if(ql==1){ facG1=0.009; facG2=0.009;   Log(2)<<"E_SymStr ";   }// sym stretch
		if(ql==2){ facG1=-0.009; facG2=-0.009; Log(2)<<"E_SymCom ";   }// sym compression	
		if(ql==3){ facG1=0.009; facG2=-0.009;  Log(2)<<"E_ASym1  ";   }// Asym 1
		if(ql==4){ facG1=-0.009; facG2=0.009;  Log(2)<<"E_ASym2  ";  }// Asym 2			
		pG1A1 -= (facG1)%bv1; // minus because bv1 points from C-O
		pG1A2 += (facG1)%bv1;
		pG2A1 -= (facG2)%bv2; // minus because bv1 points from C-O
		pG2A2 += (facG2)%bv2;
		Torder1 dG1A1A1 = pG1A1 - pG2A1;
		Torder1 dG1A1A2 = pG1A1 - pG2A2;
		Torder1 dG1A2A1 = pG1A2 - pG2A1;
		Torder1 dG1A2A2 = pG1A2 - pG2A2;
		Eexact[ql] += crd.c(g1a1)*crd.c(g2a1)/( dG1A1A1.mag() );
		Eexact[ql] += crd.c(g1a1)*crd.c(g2a2)/( dG1A1A2.mag() );
		Eexact[ql] += crd.c(g1a2)*crd.c(g2a1)/( dG1A2A1.mag() );
		Eexact[ql] += crd.c(g1a2)*crd.c(g2a2)/( dG1A2A2.mag() );
		Log(2) <<"\t"<<Eexact[ql];
		if(ql>0) Log(2)<<"\tDiff="<<Eexact[ql]-Eexact[0];
		Log(2) <<"\n";
	}
	double avg_Sym = (Eexact[1]+Eexact[2])/2 -Eexact[0];
	double avg_ASm = (Eexact[3]+Eexact[4])/2 -Eexact[0];
	Log(2) <<"Avg Sym="<<avg_Sym<<"\tAvg ASym="<<avg_ASm<<"\n";
	if(avg_Sym > avg_ASm) Log(2) <<"EXACT_INTERACTION_ORDER_LOW-HIGH_FREQ= -1,+1\n";
	if(avg_Sym < avg_ASm) Log(2) <<"EXACT_INTERACTION_ORDER_LOW-HIGH_FREQ= +1,-1\n";
	if(V12 > 0) Log(2) <<"APPROXIMATE_INTERACTION_ORDER_LOW-HIGH_FREQ= -1,+1\n";
	if(V12 < 0) Log(2) <<"APPROXIMATE_INTERACTION_ORDER_LOW-HIGH_FREQ= +1,-1\n";
	//-------------------------------------------------------------------------------------
	// calculate R
	double F = (Freq.at(nm1-1) + Freq.at(nm2-1))/2 ; // average Freq
	  Log(3)<<"cross Prod bv1 x bv2 = "<< (bv1|bv2) <<" , dot with R12 = "<< dotp(vecR12cm, (bv1|bv2) ) <<"\n";
	double R = dotp(vecR12cm, (bv1|bv2) )*D ; // R is R_plus combination
	R *= F*3.1415/2;
	R *= (-1)*1e4; // convert to (10**-44 esu**2-cm**2)
	cout <<"\n-----------------------------------------\n";
	cout <<"R+="<<R<<" F+="<<F+V12<<" cm-1\n"; 
	cout <<" R-="<<-R<<" F-="<<F-V12<<" cm-1\n";
	double dihed = xyzvec_dihedral(g1a2,g1a1,g2a1,g2a2,&crd);
	cout << "Dihedral="<<dihed<<"\n";
	double VDF_nm1 = 4*RStr.at(nm1-1)/DStr.at(nm1-1);
	double VDF_nm2 = 4*RStr.at(nm2-1)/DStr.at(nm2-1);
	double D_plus = D * (1 + dotp(bv1,bv2));
	double D_minus = D * (1 - dotp(bv1,bv2));
	double EXT_VDF_plus = 4*R/ ( D_plus );
	double EXT_VDF_minus= (-4)*R/ ( D_minus );
	
	//-------------------------------------------------------------------------------------------------------------
	Log(2)<<"TABLE_FORMAT_1:,NM_"<<nm1<<" freq,NM_"<<nm2<<" freq,NM_"<<nm1<<" DS,NM_"<<nm2<<" DS,NM_"<<nm1<<" RS,NM_"<<nm2<<" RS" \
			<<",dihedral,bv1-bv2 dot product,V12 cm-1,RS_F+,RS_F-,F+,F-\n";
	string z=", "; //delimiter
	Log(2)<<"INFO1_FORMAT_1:,"<< Freq.at(nm1-1) <<z<< Freq.at(nm2-1) <<z<< DStr.at(nm1-1) <<z<< DStr.at(nm2-1) <<z \
			<< RStr.at(nm1-1) <<z<< RStr.at(nm2-1) <<z<< dihed <<z<< DOT_PROD <<z<< V12 <<z<< R <<z<< -R <<z<< F+V12 <<z<< F-V12 <<"\n";
			
	Log(2)<<"TABLE_FORMAT_2:,NM_"<<nm1<<" freq,NM_"<<nm2<<" freq,NM_"<<nm1<<" DS,NM_"<<nm2<<" DS,NM_"<<nm1<<" RS,NM_"<<nm2<<" RS" \
			<<",dihedral,bv1-bv2 dot product,V12 cm-1,RS_F-,RS_F+,F-,F+\n";
	Log(2)<<"INFO2_FORMAT_2:,"<< Freq.at(nm1-1) <<z<< Freq.at(nm2-1) <<z<< DStr.at(nm1-1) <<z<< DStr.at(nm2-1) <<z \
			<< RStr.at(nm1-1) <<z<< RStr.at(nm2-1) <<z<< dihed <<z<< DOT_PROD <<z<< V12 <<z<< -R <<z<< R <<z<< F-V12 <<z<< F+V12 <<"\n";

	Log(2)<<"TABLE_FORMAT_3:,NM_"<<nm1<<" freq,NM_"<<nm2<<" freq,NM_"<<nm1<<" DS,NM_"<<nm2<<" DS,NM_"<<nm1<<" RS,NM_"<<nm2<<" RS" \
			<<",NM_"<<nm1<<" symmetry,NM_"<<nm2<<" symmetry,dihedral,bv1-bv2 dot product,V12 cm-1,RS Lower_freq,RS Higher_freq,Freq Lower,Freq Higher\n";	
	if(V12 > 0){ // F+V12 is larger
		Log(2)<<"INFO3_FORMAT_3:,"<< Freq.at(nm1-1) <<z<< Freq.at(nm2-1) <<z<< DStr.at(nm1-1) <<z<< DStr.at(nm2-1) <<z \
			<< RStr.at(nm1-1) <<z<< RStr.at(nm2-1) <<z<<nm1_symS<<z<<nm2_symS <<z<<dihed <<z<< DOT_PROD <<z<< V12 <<z<< \
			-R <<z<< R <<z<< F-V12 <<z<< F+V12 <<"\n";
	}
	else{ // flip R- and R+ ; F-V12 and F+V12
		Log(2)<<"INFO3_FORMAT_3:,"<< Freq.at(nm1-1) <<z<< Freq.at(nm2-1) <<z<< DStr.at(nm1-1) <<z<< DStr.at(nm2-1) <<z \
			<< RStr.at(nm1-1) <<z<< RStr.at(nm2-1) <<z<<nm1_symS<<z<<nm2_symS <<z<< dihed <<z<< DOT_PROD <<z<< V12 <<z \
			<< R <<z<< -R <<z<< F+V12 <<z<< F-V12 <<"\n";	
	}
	
	Log(2)<<"TABLE_FORMAT_4:,NM_"<<nm1<<" freq,NM_"<<nm2<<" freq,NM_"<<nm1<<" DS,NM_"<<nm2<<" DS,NM_"<<nm1<<" RS,NM_"<<nm2<<" RS" \
			<<",NM_"<<nm1<<" g,NM_"<<nm2<<" g,NM_"<<nm1<<" sym,NM_"<<nm2<<" sym,dihed,bv1-bv2 ang,V12 cm-1,RS Lower_freq,RS Higher_freq"\
			<<",Freq Lower,Freq Higher,g Low Freq, g high Freq\n";	
	if(V12 > 0){ // F+V12 is larger
		Log(2)<<"INFO4_FORMAT_4:,"<< Freq.at(nm1-1) <<z<< Freq.at(nm2-1) <<z<< DStr.at(nm1-1) <<z<< DStr.at(nm2-1) <<z \
			<< RStr.at(nm1-1) <<z<< RStr.at(nm2-1) <<z<< VDF_nm1 <<z<< VDF_nm2 <<z \
			<<nm1_symS<<z<<nm2_symS <<z<<dihed <<z<< DOT_PROD <<z<< V12 <<z \
			<< -R <<z<< R <<z<< F-V12 <<z<< F+V12 <<z<< EXT_VDF_minus <<z<< EXT_VDF_plus << "\n";
	}
	else{ // flip R- and R+ ; F-V12 and F+V12
		Log(2)<<"INFO4_FORMAT_4:,"<< Freq.at(nm1-1) <<z<< Freq.at(nm2-1) <<z<< DStr.at(nm1-1) <<z<< DStr.at(nm2-1) <<z \
			<< RStr.at(nm1-1) <<z<< RStr.at(nm2-1) <<z <<VDF_nm1 <<z<< VDF_nm2 <<z \
			<<nm1_symS<<z<<nm2_symS <<z<< dihed <<z<< DOT_PROD <<z<< V12 <<z \
			<< R <<z<< -R <<z<< F+V12 <<z<< F-V12 <<z<< EXT_VDF_plus <<z<< EXT_VDF_minus <<"\n";	
	}
	//--- geometry information output
	Torder1 R12norm = vecR12;
	R12norm.norm();
	double DOT_PROD1 = acos(dotp(bv1,R12norm))*180/3.141592;
	double DOT_PROD2 = acos(dotp(bv2,R12norm))*180/3.141592;
	Log(2)<<"TABLE_FORMAT_5:, R12 (Angstrom), R12 (cm), R12-Mu1 Ang, R12-Mu2 Ang, bv1-bv2 ang, dihed, V12 (cm-1), V12_term1 (cm-1), V12_term2 (cm-1)\n";
	Log(2)<<"INFO5_FORMAT_5:,"<< R12 <<z<< R12cm <<z<< DOT_PROD1 <<z<< DOT_PROD2 <<z<< DOT_PROD <<z<< dihed <<z \
			<< V12 <<z<< V12_t1/(29979245800*(6.6260755e-27)) <<z<< V12_t2/(29979245800*(6.6260755e-27)) <<"\n";

    Log(2)<<"TABLE_FORMAT_6:,NM_"<<nm1<<" freq,NM_"<<nm2<<" freq,NM_"<<nm1<<" DS,NM_"<<nm2<<" DS,NM_"<<nm1<<" RS,NM_"<<nm2<<" RS" \
			<<",Freq Lower,Freq Higher,DS Lower_freq,DS Higher_freq,RS Lower_freq,RS Higher_freq\n";	
	if(V12 > 0){ // F+V12 is larger
		Log(2)<<"INFO6_FORMAT_6:,"<< Freq.at(nm1-1) <<z<< Freq.at(nm2-1) <<z<< DStr.at(nm1-1) <<z<< DStr.at(nm2-1) <<z \
			<< RStr.at(nm1-1) <<z<< RStr.at(nm2-1) <<z<< F-V12 <<z<< F+V12 <<z<< D_minus <<z<< D_plus <<z<< -R <<z<< R <<"\n";	
	}
	else{ // flip R- and R+ ; F-V12 and F+V12
		Log(2)<<"INFO6_FORMAT_6:,"<< Freq.at(nm1-1) <<z<< Freq.at(nm2-1) <<z<< DStr.at(nm1-1) <<z<< DStr.at(nm2-1) <<z \
			<< RStr.at(nm1-1) <<z<< RStr.at(nm2-1) <<z<< F+V12 <<z<< F-V12 <<z<< D_plus <<z<< D_minus <<z<< R <<z<< -R <<"\n";
	}
	 
	if(Freq.at(nm1-1) > Freq.at(nm2-1) ){
			if(loudness >1 ){
			  cerr << "must specify nm1 to be lower freqency to get table 7 to print.\n";
			}
	}
	else{
		
		
	}
	cout <<"\n\
Please cite \"Determination of the Absolute Configurations Using\n\
   Exciton Chirality Method for Vibrational Circular Dichroism:\n\
   Right Answers for the Wrong Reasons?, C. L. Covington, V. P. Nicu,\n\
   P.L.Polavarapu,  J. Phys. Chem. A, 2015, 119 (42), pp 10589-10601\"\n";
			
return 0;
}
