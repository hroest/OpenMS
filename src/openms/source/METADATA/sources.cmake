### the directory name
set(directory source/METADATA)

### list all filenames of the directory here
set(sources_list
AbsoluteQuantitationStandards.cpp
DocumentIDTagger.cpp
DocumentIdentifier.cpp
ExperimentalDesign.cpp
ExperimentalSettings.cpp
MSQuantifications.cpp
PeptideEvidence.cpp
PeptideHit.cpp
PeptideIdentification.cpp
ProteinHit.cpp
ProteinIdentification.cpp
SpectrumLookup.cpp
SpectrumMetaDataLookup.cpp
SpectrumSettings.cpp
Tagging.cpp
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
Acquisition.cpp
AcquisitionInfo.cpp
ContactPerson.cpp
DataArrays.cpp
Digestion.cpp

ChromatogramSettings.cpp
Precursor.cpp

Identification.cpp
IdentificationHit.cpp
Instrument.cpp
InstrumentSettings.cpp

Modification.cpp

CVTerm.cpp
CVTermList.cpp
CVTermListInterface.cpp
MetaInfo.cpp
MetaInfoDescription.cpp
MetaInfoInterface.cpp
MetaInfoRegistry.cpp
Software.cpp
DataProcessing.cpp
Sample.cpp
SampleTreatment.cpp
ScanWindow.cpp
SourceFile.cpp
IonDetector.cpp
IonSource.cpp
MassAnalyzer.cpp
Gradient.cpp
HPLC.cpp
Product.cpp

SpectrumIdentification.cpp
)

### add path to the filenames
set(sources)
foreach(i ${sources_list})
	list(APPEND sources ${directory}/${i})
endforeach(i)

### pass source file list to the upper instance
set(OpenMS_BASE_sources ${OpenMS_BASE_sources} ${sources})

### source group definition
source_group("Source Files\\METADATA" FILES ${sources})

