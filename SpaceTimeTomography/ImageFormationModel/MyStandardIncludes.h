#ifndef __MyStandardIncludes_h
#define __MyStandardIncludes_h

#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__NT__)
	#define My_OS_WIN32
#elif defined(__linux__) || defined(__linux)
	#define My_OS_LINUX
#endif

#if defined(My_OS_WIN32) 
#define _WIN32_WINNT 0x0501 
#include <windows.h>
	#if defined(My_STANDARD_DLL)  
		#define My_STANDARD_API __declspec(dllexport)
	#else
		#define My_STANDARD_API __declspec(dllimport)
	#endif
#else   //Other system
#undef My_STANDARD_DLL
#endif

#ifndef My_STANDARD_API
	#define My_STANDARD_API
#endif

#endif 