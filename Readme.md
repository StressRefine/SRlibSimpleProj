SRlibSimpleProj

Updated Version will compile on linux or windows
linux Version: download source, go into folder linuxDebug or linuxRelease and make all there, library will reside thrre
note that both the makefile and all ".mk" files are needed because the makefile references them

Project to create the stressRefine library SRlibSimple
Windows Visual Studio version that was last modified using VS2013.
That and any newer version should work to compile this project.
It is in C++
Doubleclick on the ".sln" file to open the project.
If using a newer version than 2013, you may get a message about updating to the latest VS version, to which you should say yes.

Choose a configuration at the top (debug or release). Next to that it should say x64.
You can change that to win32 but executation is faster in x64 mode.
Now hit build to compile. A subfolder is created x64/debug or x64/release depending on the configuration.
Inside that folder will be the library file SRlibSimple.lib
