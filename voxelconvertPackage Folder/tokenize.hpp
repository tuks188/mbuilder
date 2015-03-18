#ifndef TOKENIZE_HPP
#define TOKENIZE_HPP

#include<vector>
#include<string>
using namespace std;
/*
  void tokenize(const string& str, std::vector<string>& tokens, const
  string& delimiters = " ")
  ====================================================================
  Taken from http://www.oopweb.com/CPP/Documents/CPPHOWTO/Volume/C++Programming-HOWTO-7.html

  "A very common operation with strings, is to tokenize it with a
  delimiter of your own choice. This way you can easily split the
  string up in smaller pieces, without fiddling with the find()
  methods too much. In C, you could use strtok() for character arrays,
  but no equal function exists for strings. This means you have to
  make your own. Here is a couple of suggestions, use what suits your
  best."
  ====================================================================
*/
void tokenize(const string& str,
                      std::vector<string>& tokens,
                      const string& delimiters = " ")
{
    // Skip delimiters at beginning.
    string::size_type lastPos = str.find_first_not_of(delimiters, 0);

    // Find first "non-delimiter".
    string::size_type pos     = str.find_first_of(delimiters, lastPos);

    while (string::npos != pos || string::npos != lastPos)
    {
        // Found a token, add it to the vector.
        tokens.push_back(str.substr(lastPos, pos - lastPos));

        // Skip delimiters.  Note the "not_of"
        lastPos = str.find_first_not_of(delimiters, pos);

        // Find next "non-delimiter"
        pos = str.find_first_of(delimiters, lastPos);
    }
}




#endif
