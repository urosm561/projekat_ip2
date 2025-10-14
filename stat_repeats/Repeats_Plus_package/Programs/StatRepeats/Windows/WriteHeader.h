#ifndef WRITEHEADER_H 
#define WRITEHEADER_H 


#include <iostream>
#include "CommandLine.h"

void WriteHeader(const CommandLine& c, FILE* writer);

void WriteHeader(const CommandLine& c, std::ostream& writer);
#endif