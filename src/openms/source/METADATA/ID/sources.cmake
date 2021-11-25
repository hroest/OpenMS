### the directory name
set(directory source/METADATA/ID)

### list all filenames of the directory here
set(sources_list
# IdentificationData.cpp
IdentificationDataConverter.cpp
)

### add path to the filenames
set(sources)
foreach(i ${sources_list})
	list(APPEND sources ${directory}/${i})
endforeach(i)

### pass source file list to the upper instance
set(OpenMS_CORE_sources ${OpenMS_CORE_sources} ${sources})
set(OpenMS_BASE_sources ${OpenMS_BASE_sources} ${directory}/IdentificationData.cpp)

### source group definition
source_group("Source Files\\METADATA\\ID" FILES ${sources})
