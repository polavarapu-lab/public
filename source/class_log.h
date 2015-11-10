
#ifndef CLASS_LOG_
#define CLASS_LOG_

extern int loudness; //set with -l  0 very quite, 1 normal, 2 noisy, 3 debug , -1 error
class Log //custom print message class. .Use like   Log()<<"message \n"; or   Log(loudness,FN) << "message \n";
{
	int loud;
public:
	Log(){
		loud =1;
	}
    Log(int a){	
        loud = a;
    }
	Log(int a, const std::string &funcName){	
        loud = a;
		if(a>2) std::cout << funcName <<": ";
    }
	Log(int a, const std::string &funcName, int err){
		loud =a;  //for common errortypes see below
		if(a==-2) std::cerr <<"Warning Type "<<err<<" Thrown in "<<funcName<<": ";
		if(a==-1) std::cerr <<"Error Type "<<err<<" Thrown in "<<funcName<<": ";  
	}
    template <class T>
    Log &operator<<(const T &v)
    {
		if(loud == -1){
			std::cerr << v ;
		}
		else if(loud <= loudness){
			std::cout << v ;
		}
        return *this;
    }

    ~Log(){}
};

#endif

/* common error types
	0= catch exception
	1=file load
	2=math warning
	3=user input error
*/


/* above was updated on 7-30-2014 from old which had weird error output

old one

extern int loudness; //set with -l  0 very quite, 1 normal, 2 noisy, 3 debug , -1 error
class Log //custom print message class. .Use like   Log()<<"message \n"; or   Log(loudness,FN) << "message \n";
{
        int loud;
public:
        Log(){
                loud =1;
        }
    Log(int a){
        loud = a;
    }
        Log(int a, const std::string &funcName){
        loud = a;
                if(a>2) std::cout << funcName <<": ";
    }
    template <class T>
    Log &operator<<(const T &v)
    {
                if(loud == -1){
                        std::cerr << "ERROR:" << v << endl;
                }
                if(loud <= loudness){
                        std::cout << v ;
                }
        return *this;
    }

    ~Log(){}
};




*/