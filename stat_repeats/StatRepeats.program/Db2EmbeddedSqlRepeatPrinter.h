#ifndef DB2EMBEDDEDSQLREPEATPRINTER_H
#define DB2EMBEDDEDSQLREPEATPRINTER_H

#ifdef DB2EMBEDDEDOUTPUT

#include <string>
#include "RepeatPrinter.h"
#include "CommandLine.h"

RepeatPrinter* CreateDb2EmbeddedSqlRepeatPrinter(const CommandLine& c);

#endif

#endif