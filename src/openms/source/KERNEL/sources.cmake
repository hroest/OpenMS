### the directory name
set(directory source/KERNEL)

### list all filenames of the directory here
set(sources_list
AreaIterator.cpp
ConversionHelper.cpp
MRMFeature.cpp
MRMTransitionGroup.cpp
OnDiscMSExperiment.cpp
PeakIndex.cpp
RangeManager.cpp
StandardTypes.cpp
ChromatogramTools.cpp
SpectrumHelper.cpp
)

### add path to the filenames
set(sources)
foreach(i ${sources_list})
	list(APPEND sources ${directory}/${i})
endforeach(i)

### pass source file list to the upper instance
set(OpenMS_CORE_sources ${OpenMS_CORE_sources} ${sources})

### list all filenames of the directory here
set(sources_list
MSExperiment.cpp
MSChromatogram.cpp
MSSpectrum.cpp
MassTrace.cpp
RichPeak2D.cpp
ChromatogramPeak.cpp
DPeak.cpp
Peak1D.cpp
Peak2D.cpp
BaseFeature.cpp # PeptideIdentification
ConsensusFeature.cpp
ConsensusMap.cpp
Feature.cpp
FeatureHandle.cpp
FeatureMap.cpp
)

### add path to the filenames
set(sources)
foreach(i ${sources_list})
	list(APPEND sources ${directory}/${i})
endforeach(i)

### pass source file list to the upper instance
set(OpenMS_BASE_sources ${OpenMS_BASE_sources} ${sources})

### source group definition
source_group("Source Files\\KERNEL" FILES ${sources})

