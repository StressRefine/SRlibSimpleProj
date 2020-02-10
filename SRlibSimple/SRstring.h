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
// SRstring.h: interface for the SRstring class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(SRSTRING_INCLUDED)
#define SRSTRING_INCLUDED

#include <string.h>

static char bdfBuf[17];

class SRstring
{
public:
	void Left(char c, SRstring &s2, bool last = true);
	void Left(int n, SRstring &s2);
	void Right(char c, SRstring &s2);
	void Right(int n, SRstring &s2);
	void Copy(char* s, int n = 0);
	void Copy(SRstring& s2){ Copy(s2.str); };
	void Cat(char* s, int n = 0);
	void Cat(SRstring& s2){ Cat(s2.str); };
	void operator = (char* s) { Copy(s); };
	void operator = (SRstring& s2) { Copy(s2.str); };
	void operator = (SRstring* s2) { Copy(s2->str); };
	void operator += (char* s) { Cat(s); };
	void operator += (SRstring& s) { Cat(s); };
	bool Compare(char* s2, int n = 0);
	bool Compare(SRstring& s2, int n = 0){ return Compare(s2.str, n); };
	bool CompareCaseSensitive(char* s2, int n = 0);
	bool CompareUseLength(char* s2, bool useCase = false);
	bool SRstring::CompareSkipBlanks(char *s2);
	bool CompareUseLength(SRstring& s2, bool useCase = false) { return CompareUseLength(s2.str, useCase); };
	bool isCommentOrBlank(bool skipContinuation = false);
	bool operator == (char* s2) { return CompareUseLength(s2); };
	bool operator != (char* s2) { return !Compare(s2); };
	bool operator == (SRstring& s2) { return CompareUseLength(s2); };
	char* Token(char* sep = NULL);
	bool TokRead(int& i);
	bool TokRead(double& r,  bool checkForTrailingComment = false);
	double RealRead();
	int IntRead();
	void setTokSep(char *sep){ tokSep = sep; };
	bool continueCheck();

	char GetChar(int i){ return str[i]; };
	char operator [] (int i) { return GetChar(i); };
	//FirstChar finds 1st occurrence of character c. returns c and remainder of
	//string to right of c; returns NULL if c not found
	char* FirstChar(char c) { return strchr(str, c); };
	int FirstCharLocation(char c);
	void truncate(int n);
	//LastChar finds last occurrence of character c. returns c and remainder of
	//string to right of c; returns NULL if c not found
	char* LastChar(char c, bool after = false)
	{
		if (after)
			return (strrchr(str, c) + 1);
		else
			return strrchr(str, c);
	};
	void Clear();
	bool isAllBlank();

	SRstring(SRstring& s2) { fresh = true; len = strlen(s2.str); str = new char[len + 1]; strcpy_s(str, len + 1, s2.str); nextToken = NULL; tokSep = NULL; };
	SRstring(char* s) { fresh = true; len = strlen(s); str = new char[len + 1]; strcpy_s(str, len + 1, s); nextToken = NULL; tokSep = NULL; };
	SRstring() { fresh = true; len = 0; str = (char *)0; tokSep = NULL; csv = false; };
	~SRstring() { Clear(); };

	int len;
	char* str;
	int tokindex;
	char *nextToken;
	char *tokSep;
	bool csv;
protected:
	bool fresh;
};

#endif //if !defined(SRSTRING_INCLUDED)
