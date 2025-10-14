#ifndef COMMAND_LINE_PARSER_H 
#define COMMAND_LINE_PARSER_H 

#include <string>
#include <vector>
#include <unordered_map>
#include <tuple>
#include <stdexcept>


std::unordered_map<std::string, std::string> ParseCommandLine(int argsLength, char** args, const std::initializer_list<std::string>& mandatoryArguments, const std::initializer_list<std::string>& otherArguments, std::initializer_list<std::tuple<std::string, std::string, std::string>>& shortcuts);
bool CaseInsensitiveStartsWith(const std::string& str1, const std::string& str2); 

#endif