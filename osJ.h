//To get Windows to export symbols from .dll library

//To get Windows to export symbols from .dll library
#ifdef _WIN32
#define EXP __declspec(dllexport)
#else
#define EXP
#include <stdio.h>
inline bool fopen_s(FILE** file,const char* fname,const char* rw){
	*file = fopen(fname, rw);
	if(*file) return true;
	else return false;
}
#endif
