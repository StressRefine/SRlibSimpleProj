/*
Copyright (c) 2020 Richard King

The StressRefine library is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

The StressRefine library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

The terms of the GNU General Public License are explained in the file COPYING,
also available at <https://www.gnu.org/licenses/>

*/

//////////////////////////////////////////////////////////////////////
//
// SRstring.cpp: implementation of the SRstring class.
//
//////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <string.h>
#include "SRstring.h"

void SRstring::Copy(char *s, int n)
{
//Copy string s, overriding SRstring.str of "this"
	//input:
		//s = string
		//n = number of characters to copy. If n is 0, copy length of s

	Clear();
	if(s == NULL)
		return;
	if(n == 0)
		len = strlen(s);
	else
		len = n;
	str = new char[len + 1];
	if(n == 0)
		strcpy_s(str, len + 1, s);
	else
	{
		strncpy_s(str, len + 1, s, n);
		str[len] = '\0';
	}
	fresh = true;
}

void SRstring::Cat(char *s, int n)
{
//concatenate string s onto SRstring.str, reallocating space
	//input:
		//s = string
		//n = number of characters to concat. If n is 0, Copy s

	if(len == 0)
	{
		Copy(s);
		return;
	}
	else
	{
		char *t = new char[len + 1];
		strcpy_s(t, len + 1, str);
		delete str;
		if(n == 0)
			len += strlen(s);
		else
			len += n;
		str = new char[len + 1];
		strcpy_s(str, len + 1, t);
		if(n == 0)
			strcat_s(str, len + 1, s);
		else
			strncat_s(str, len + 1, s, n);
		delete t;
	}
	fresh = true;
}

bool SRstring::Compare(char *s2, int n)
{
//see if s2 is same as str ("stricmp"); optional n chars ("strnicmp")
//ignoring case

	bool ret = false;
	if(len == 0)
	{
		if (strcmp(s2, "") == 0)
			return true;
		else
			return false;
	}
	if(n == 0)
	{
		if(_stricmp(str,s2) == 0)
			ret = true;
	}
	else
	{
		if(_strnicmp(str,s2,n) == 0)
			ret = true;
	}
	return ret;
}

bool SRstring::CompareCaseSensitive(char *s2, int n)
{
//see if s2 is same as str ("strcmp"); optional n chars ("strncmp")
//case sensitive

	bool ret = false;
	if(len == 0)
	{
		if(strcmp(s2,"") == 0)
			return true;
		else
			return false;
	}
	if(n == 0)
	{
		int comp = strcmp(str, s2);
		if (comp == 0)
			ret = true;
	}
	else
	{
		if(strncmp(str,s2,n) == 0)
			ret = true;
	}
	return ret;
}

char *SRstring::Token(char *sep)
{
//find next "Token" in SRstring.str, using standard library function strtok
//caution: like strtok, this destroys the original string as it works its way through it 
//finding tokens

	if (len == 0)
		return false;
	char *mysep;
	if (sep == NULL)
	{
		if (tokSep != NULL)
			mysep = tokSep;
		else
			mysep = " ";
	}
	else
		mysep = sep;
	if(fresh)
	{
		fresh = false;
		return strtok_s(str, mysep, &nextToken);
	}
	else
	{
		return strtok_s(NULL, mysep, &nextToken);
	}
}

bool SRstring::TokRead(int &i)
{
//get next Token of this string, interpret as integer

	char *s = Token();
	if (s == NULL)
		return false;
	if (sscanf_s(s, "%d", &i) > 0)
		return true;
	else
		return false;
}

bool SRstring::TokRead(double  &r, bool checkForTrailingComment)
{
//get next Token of this string, interpret as double

	char *s = Token();
	r = 0.0;
	if (s == NULL)
		return false;
	if (checkForTrailingComment)
	{
		SRstring tmp = s;
		if (tmp.Compare("//", 2))
			return false;
	}
	if (sscanf_s(s, "%lg", &r) > 0)
		return true;
	else
		return false;
}

int SRstring::IntRead()
{
	//read an integer from this string

	int i = 0;
	sscanf_s(str, "%d", &i);
	return i;
}

double SRstring::RealRead()
{
	//read a double from this string

	double r = 0.0;
	sscanf_s(str, "%lg", &r);
	return r;
}

void SRstring::Clear()
{
	//clear a string for reuse

	if (len != 0)
	{
		delete str;
		str = NULL;
	}
	len = 0;
	fresh = true;
};

bool SRstring::CompareSkipBlanks(char *s2)
{
	//see if s2 is same as str ("strncmp") using length of s2
	//as number of characters

	int n = strlen(str);
	int nskip = 0;
	for (int i = 0; i < n; i++)
	{
		if (str[i] != ' ')
			break;
		nskip++;
	}
	char *s3 = str + nskip;
	n = strlen(s2);
	if (_strnicmp(s2, s3, n) == 0)
		return true;
	else
		return false;
}

bool SRstring::CompareUseLength(char *s2, bool useCase)
{
//see if s2 is same as str ("strncmp") using length of s2
//as number of characters

	if(len == 0)
	{
		if(strcmp(s2,"") == 0)
			return true;
		else
			return false;
	}
	int n = strlen(s2);
	if(n == 0 && len != 0)
		return false;
	if(useCase)
	{
		if (strncmp(str, s2, n) == 0)
			return true;
		else
			return false;
	}
	else
	{
		if (_strnicmp(str, s2, n) == 0)
			return true;
		else
			return false;
	}
}

void SRstring::Left(int n, SRstring &s2)
{
	//copy leftmost n characters of "this" to s2

	s2.Copy(str, n);
}

void SRstring::Left(char c, SRstring &s2, bool last)
{
	//copy characters of "this" left of char c to s2
	//if last = true, then use last occurrence of c, else use 1st

	int n;
	char *tmp;
	if(last)
		tmp = LastChar(c);
	else
		tmp = FirstChar(c);
	if (tmp != NULL)
	{
		//integer location in string:
		n = tmp - str + 1;
		Left(n - 1, s2);
	}
}

void SRstring::Right(char c, SRstring &s2)
{
	//copy characters of "this" right of last occurrence of char c to s2

	int n;
	char *tmp;
	tmp = LastChar(c);
	if (tmp != NULL)
	{
		//integer location in string:
		n = tmp - str + 1;
		s2.Copy(str + n);
	}
}

void SRstring::Right(int n, SRstring &s2)
{
	//copy rightmost chars of this string to S2

	s2.Copy(str + n);
}

bool SRstring::isCommentOrBlank(bool skipContinuation)
{
	//check if this string is a comment line or blank
	//input:
		//skipContinuation = true to also skip continuation character ("#") in multi-line input

	if (Compare("") || CompareUseLength("\\") || CompareUseLength("//") || (skipContinuation && CompareUseLength("#")))
		return true;
	else
		return false;

}

void SRstring::truncate(int n)
{
	//truncate this field at char n
	str[n] = '\0';
	len = n;
}

int SRstring::FirstCharLocation(char c)
{
	//returns 1st location of char c, -1 if not found
	char *cloc;
	int loc = -1;

	// Search forward. 
	cloc = strchr(str, c);
	if (cloc != NULL)
		loc = (int)(cloc - str + 1);
	return loc;
}

bool SRstring::isAllBlank()
{
	for (int i = 0; i < len; i++)
	{
		if (str[i] != ' ')
			return false;
	}
	return true;
}


