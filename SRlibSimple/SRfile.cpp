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
// SRfile.cpp: implementation of the SRfile class
//
//////////////////////////////////////////////////////////////////////

#include <stdarg.h>
#include <stdio.h>
#include <string.h>
#include <io.h>
#include <direct.h>
#include "SRfile.h"
#include "SRmodel.h"
#include "SRmachDep.h"

#define ECHOTOSTAT true

extern SRmodel model;

bool SRfile::Existcheck(char *name)
{
    //determine if file "name" exists
    //return:
		//true if it exists else false

	if(_access(name,0) != -1)
		return true;
	else
		return false;
}

bool SRfile::Existcheck()
{
	//determine if a file corresponding to class variable filename exists

	return SRfile::Existcheck(filename.str);
}

bool SRfile::Open(FileOpenMode mode, char *name)
{
    //open file "name"
    //input:
        //mode = SRinputMode, SRoutputMode, SRappendMode, SRoutbinaryMode, SRinbinaryMode, or SRinoutbinaryMode
    //return:
		//true if file was already opened else false
	if(opened)
		return false;
	filePos = 0; //position for binary i/o is at beginning of file for newly opened file

	if (name != NULL)
		filename = name;
	if (filename.len == 0)
		return false;

	if (mode == SRinputMode)
	{
		//unsuccessful trying to open a non-existent file for reading
		if (!Existcheck(filename.str))
			return false;
		fopen_s(&fileptr, filename.str, "r");
	}
	else if (mode == SRoutputMode)
	{
		//overwrite existing file for output mode:
		if (Existcheck(filename.str))
			_unlink(filename.str);
		fopen_s(&fileptr, filename.str, "w");
	}
	else if (mode == SRappendMode)
		fopen_s(&fileptr, filename.str, "a");
	else if (mode == SRoutbinaryMode)
		fopen_s(&fileptr, filename.str, "wb");
	else if (mode == SRinbinaryMode)
	{
		//unsuccessful trying to open a non-existent file for reading
		if (!Existcheck(filename.str))
			return false;
		fopen_s(&fileptr, filename.str, "rb");
	}
	else if (mode == SRinoutbinaryMode)
		fopen_s(&fileptr, filename.str, "rb+");
	else
		return false;

	if (fileptr == NULL)
	{
		SRASSERT(0);
		return false;
	}
	else
	{
		opened = true;
		return true;
	}
}

bool SRfile::GetLine(SRstring &line,bool noSlashN)
{
    //get a line from a file
    //input:
        //noSlashN = true to not return the "\n" character at end of line else false
    //output:
        //line = the fetched line stored as SRstring
    //return
        //true if successful else false (e.g. EOF)
	if (!opened)
		return false;
	line.Clear();
	char *tmp,c;
	int len;
	tmp = fgets(linebuf, MAXLINELENGTH, fileptr);
	if (tmp == NULL)
		return false;
	if(noSlashN)
	{
		len = strlen(tmp);;
		c = tmp[len-1];
		if(c == '\n')
			tmp[len-1] = '\0';
	}
	else
	{
		len = strlen(tmp);
		c = tmp[len-1];
		if(c != '\n')
		{
			//append \n if not already there:
			tmp[len] = '\n';
		}
	}
	line = tmp;
	return true;
}

bool SRfile::Close()
{
    //close a file
    //return:
		//false if file was not opened or close is unsuccessful, else true

	if(!opened)
		return false;
	opened = false;
	if(fclose(fileptr) != 0)
		return false;
	fileptr = NULL;
	return true;
}

bool SRfile::Print(char *fmt,...)
{
	//print to file; layer over vfprintf
	//input:
		//fmt = format string
		//... = variable input

	va_list arglist;
	va_start(arglist, fmt);
	int ret = vfprintf(fileptr, fmt, arglist);
	va_end(arglist);
	if(ret < 0)
		return false;
	else
		return true;
}

bool SRfile::VPrint(char *fmt, va_list arglist)
{
	//print to file; layer over vfprintf
	//input:
		//fmt = format string
		//va_list = variable argument list

	int len = strlen(fmt);
	int ret = 0;
	ret = vfprintf(fileptr, fmt, arglist);
	va_end(arglist);
	if (ret < 0)
		return false;
	else
		return true;
}

bool SRfile::PrintLine(char *fmt,...)
{
	//print to file; layer over vfprintf. append \n
	//input:
		//fmt = format string
		//... = variable input

	va_list arglist;
	va_start(arglist, fmt);
	bool ret = VPrintLine(fmt, arglist);
	va_end(arglist);
	return ret;
}

bool SRfile::VPrintLine(char *fmt, va_list arglist)
{
	//print to file; layer over vfprintf. append \n
	//input:
		//fmt = format string
		//va_list = variable argument list

	int len = strlen(fmt);
	int ret = 0;
	if(fmt[len-1] != '\n')
	{
		strcpy_s(linebuf, MAXLINELENGTH, fmt);
		strcat_s(linebuf, MAXLINELENGTH, "\n");
		ret = vfprintf(fileptr, linebuf, arglist);
	}
	else
		ret = vfprintf(fileptr, fmt, arglist);
	va_end(arglist);
	if(ret < 0)
		return false;
	else
		return true;
}

bool SRfile::PrintReturn()
{
    //print "\n" to a file
	int ret = fprintf(fileptr,"\n");
	if(ret < 0)
		return false;
	else
		return true;
}

//static:
void SRfile::GetCurrentDir(SRstring &dir)
{
    //get the current working directory for io
	//output:
		//dir = string, name of current working directory
	char buf[MAXLINELENGTH];
	char *s = _getcwd(buf, MAXLINELENGTH);
	dir = s;
}

void SRfile::Delete(char *name)
{
    //delete file "name"
	if(Existcheck(name))
	    _unlink(name);
}

void SRfile::Delete()
{
    //delete this file
	if(opened)
		Close();
	SRfile::Delete(filename.str);
}

void SRfile::PrintReportFile(char *fmt, ...)
{
	//print to model reportFile and screen using format fmt; append \n
	//input:
		//fmt = format string
		//... variable data to print
	//return:
		//true if successful else false

	SRfile* f = &model.reportFile;
	f->Open(SRappendMode);
	va_list arglist;
	va_start(arglist, fmt);
	bool ret = f->VPrintLine(fmt, arglist);
	f->Close();
}

void SRfile::PrintReportFileNoRet(char *fmt, ...)
{
	//print to model reportFile and screen using format fmt; do not append \n
	//input:
		//fmt = format string
		//... variable data to print
	//return:
		//true if successful else false

	SRfile* f = &model.reportFile;
	f->Open(SRappendMode);
	va_list arglist;
	va_start(arglist, fmt);
	bool ret = f->VPrint(fmt, arglist);
	f->Close();
}

void SRfile::LogPrint(char *fmt, ...)
{
	//print to model logFile and screen using format fmt; append \n
	//input:
		//fmt = format string
		//... variable data to print
	//return:
		//true if successful else false

	SRfile* f = &model.logFile;
	f->Open(SRappendMode);
	va_list arglist;
	va_start(arglist, fmt);
	bool ret = f->VPrintLine(fmt, arglist);
	f->Close();
}

void SRfile::LogPrintNoRet(char *fmt, ...)
{
	//print to model logFile and screen using format fmt; do not append \n
	//input:
		//fmt = format string
		//... variable data to print
	//return:
		//true if successful else false

	SRfile* f = &model.logFile;
	f->Open(SRappendMode);
	va_list arglist;
	va_start(arglist, fmt);
	bool ret = f->VPrint(fmt, arglist);
	f->Close();
}

void SRfile::OutPrint(char *fmt, ...)
{
	//print to model output File using format fmt; append \n
	//input:
		//fmt = format string
		//... variable data to print
	//return:
		//true if successful else false

	SRfile* f = &model.outFile;
	f->Open(SRappendMode);
	va_list arglist;
	va_start(arglist, fmt);
	bool ret = f->VPrintLine(fmt, arglist);
	f->Close();
}

void SRfile::OutPrintNoRet(char *fmt, ...)
{
	//print to model output File using format fmt; do not append \n
	//input:
		//fmt = format string
		//... variable data to print
	//return:
		//true if successful else false

	SRfile* f = &model.outFile;
	f->Open(SRappendMode);
	va_list arglist;
	va_start(arglist, fmt);
	bool ret = f->VPrint(fmt, arglist);
	f->Close();
}

bool SRfile::OutOpenNoFail()
{
	//Open a file whose name was previously assigned for output
	//return true if open successful else false
	//overwrite existing file for output mode:

	if (Existcheck(filename.str))
		_unlink(filename.str);
	fopen_s(&fileptr, filename.str, "w");
	if (fileptr == NULL)
		return false;
	else
	{
		opened = true;
		return true;
	}
}

bool SRfile::CreateDir(char *name)
{
    //create directory "name"
    //return:
		//true if successful else false

	bool ret;
	if(Existcheck(name))
		ret = true;
	else
		ret = SRmachDep::CreateDir(name);
	SRASSERT(ret);
	return ret;
}

bool SRfile::Rename(SRstring& newName, bool overWriteOk)
{
	SRstring oldname = filename;
	if (overWriteOk)
		Delete(newName.str);
	filename = newName;
	int status = rename(oldname.str, filename.str);
	if (status == 0)
		return true;
	else
		return false;
}

bool SRfile::Rename(char* newName, bool overWriteOk)
{
	SRstring s;
	s = newName;
	return Rename(s, overWriteOk);
}

bool SRfile::Screenprint(char *fmt, ...)
{
	//print to cmd screen
	//input:
		//fmt = format string
		//... variable data to print
	//return:
		//true if successful else false
	va_list arglist;
	va_start(arglist, fmt);
	return vprintf_s(fmt, arglist);
}


