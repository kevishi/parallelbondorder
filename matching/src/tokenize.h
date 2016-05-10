// tokenizing method

#ifndef TOKENIZE_H
#define TOKENIZE_H

#include <string>
#include <vector>

  //! Break a string into tokens
  //! Tokens are determined by the delimiters supplied
  //! (defaults to whitespace (i.e., spaces, tabs, newlines)
std::vector<std::string> tokenize(const std::string buffer, const char *delim = " \t\n\r");

#endif
