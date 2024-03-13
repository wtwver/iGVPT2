#include "Utils.h"

#include <sstream>
#include <algorithm>
#include <iterator>

#define debug  false
/************************************************************************************************************************************/
vector<string> split(string str)
{
        istringstream iss(str);
        vector<string> tokens;
        copy(istream_iterator<string>(iss), istream_iterator<string>(), back_inserter<vector<string> >(tokens));
        return tokens;
}
/************************************************************************************************************************************/
char *getSuffixNameFile(const char* allname)
{
	char* name = strdup(allname);
	int len=strlen(allname);
	int i;
	for(i=len;i>0;i--)
	if(name[i]=='.')
	{
		name[i] = '\0';
		break;
	}
	return name;
}
