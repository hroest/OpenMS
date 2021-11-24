### the directory name
set(directory source/SYSTEM)

### list all filenames of the directory here
set(sources_list
ExternalProcess.cpp
FileWatcher.cpp
UpdateCheck.cpp
NetworkGetRequest.cpp
PythonInfo.cpp
JavaInfo.cpp
StopWatch.cpp
SysInfo.cpp
RWrapper.cpp
File.cpp
)

### add path to the filenames
set(sources)
foreach(i ${sources_list})
	list(APPEND sources ${directory}/${i})
endforeach(i)

### pass source file list to the upper instance
set(OpenMS_BASE_sources ${OpenMS_BASE_sources} ${sources})

### source group definition
source_group("Source Files\\SYSTEM" FILES ${sources})

