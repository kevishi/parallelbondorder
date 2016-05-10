// tokenizing method

#include "tokenize.h"

std::vector<std::string> tokenize(const std::string buffer, const char *delim)
  {
    std::vector<std::string> results;
    if (buffer.size() == 0)
      return results; // empty

    std::string s(buffer);
    s += delim[0]; //forces last token to be parsed
    size_t start=0,end=0;

    for (;;)
      {
        start = s.find_first_not_of(delim,start);
        end = s.find_first_of(delim,start);

        if (end <= s.size() && start <= s.size())
          results.push_back(s.substr(start, end-start));
        else
          break;

        start = end+1;
      }

    return(results);
  }
