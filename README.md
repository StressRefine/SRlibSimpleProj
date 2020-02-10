# bdfTranslateProj
 Nastran inpurt file translator for stressRefine
Windows Visual Studio version
that was last modified using
VS2013. That and any newer version should work to compile them.
They are all in C++ except the User interface, srui, which is in C#.
Go in a folder and doubleclick on the ".sln" file to open the project.
If using a newer version than 2013, you may get a message about updating to the latest
VS version,
to which you should say yes.
 
Choose a configuration at the top (debug or release). Next to that it should say x64.
You can change that to win32 but executation is faster in x64 mode.
Now hit build to compile.
A subfolder is created x64/debug or x64/release depending on the configuration.
Inside that folder will the executable bdfTranslate.exe
