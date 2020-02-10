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
// SRfile.h: interface for the SRfile class.
//
//////////////////////////////////////////////////////////////////////

#if !(defined SRFILE_INCLUDED)
#define SRFILE_INCLUDED

#include "SRstring.h"
#include "SRutil.h"

#define MAXLINELENGTH 256

#define LOGPRINT SRfile::LogPrint
#define LOGPRINTNORET SRfile::LogPrintNoRet
#define LOGPRINTRET SRfile::LogPrint("\n");
#define OUTPRINT SRfile::OutPrint
#define OUTPRINTNORET SRfile::OutPrintNoRet
#define OUTPRINTRET SRfile::OutPrint("\n");
#define SCREENPRINT SRfile::Screenprint
#define REPPRINT SRfile::PrintReportFile
#define REPPRINTNORET SRfile::PrintReportFileNoRet

//log file is intended for use as a detailed log
//report file is intended for use as a terse summary


enum FileOpenMode{ SRinputMode, SRoutputMode, SRappendMode, SRoutbinaryMode, SRinbinaryMode, SRinoutbinaryMode };

class SRfile
{
public:
	bool OutOpenNoFail();
	static bool CreateDir(char* name);
	static void GetCurrentDir(SRstring& dir);
	static void Delete(char* name);
	static bool Existcheck(char* name);
	static bool Existcheck(SRstring& name){ return Existcheck(name.str); };
	static void PrintReportFile(char *fmt, ...);
	static void PrintReportFileNoRet(char *fmt, ...);
	static void LogPrint(char *fmt, ...);
	static void LogPrintNoRet(char *fmt, ...);
	static void OutPrint(char *fmt, ...);
	static void OutPrintNoRet(char *fmt, ...);
	bool Existcheck();
	bool VPrintLine(char* fmt, va_list arglist);
	bool VPrint(char* fmt, va_list arglist);
	void Delete();
	bool PrintReturn();
	bool Close();
	void ToTop(){ rewind(fileptr); };
	bool GetLine(SRstring& line, bool noSlashN = true);
	bool Open(FileOpenMode mode, char* name = NULL);
	bool Open(SRstring& fn, FileOpenMode mode){ return Open(mode, fn.str); };
	bool Print(char* s, ...);
	bool PrintLine(char* s, ...);
	void setFileName(SRstring& name){ filename = name; };
	void catFileName(SRstring& name){ filename += name; };
	void setFileName(char* name){ filename = name; };
	void catFileName(char* name){ filename += name; };
	int getFilePos(){ return filePos; };
	SRstring& getFileName(){ return filename; };
	bool Rename(SRstring& newName, bool overWriteOk = false);
	bool Rename(char* newName, bool overWriteOk = false);

	static bool SRfile::Screenprint(char *fmt, ...);

	SRfile(){ fileptr = NULL; opened = false; filename = ""; bdfLineSaved = false; };
	~SRfile(){Close();};

private:

	SRstring tmpstr;
	FILE* fileptr;
	bool opened;
	SRstring filename;
	char linebuf[MAXLINELENGTH];
	int filePos;
	bool bdfLineSaved;
	SRstring bdfLineSave;
};
#endif //if !(defined SRFILE_INCLUDED)