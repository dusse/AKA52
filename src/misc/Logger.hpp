//
//  Logger.hpp

#ifndef Logger_hpp
#define Logger_hpp

#include <stdio.h>
#include <mpi.h>
#include <set>

#define DEBUG    2
#define INFO     1
#define CRITICAL 0


#ifdef LOG
    #define MINIMAL_LEVEL 2
#else
    #define MINIMAL_LEVEL 1
#endif

class Logger
{
public:
    Logger();
    void writeMsg(const char* input);
    void writeMsg(const char* input, int level);
    ~Logger();
};
#endif
