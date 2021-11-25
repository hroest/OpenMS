### the directory name
set(directory source/METADATA)

### list all filenames of the directory here
set(sources_list
AbsoluteQuantitationStandards.cpp # FeatureMap
ExperimentalDesign.cpp # FeatureMap
MSQuantifications.cpp # FeatureMap
SpectrumLookup.cpp # MSSpectrum
SpectrumMetaDataLookup.cpp # FileHandler

ProteinHit.cpp # ResidueModification
PeptideHit.cpp # AASequence

ExperimentalSettings.cpp # ProteinIdentification
SpectrumSettings.cpp # PeptideIdentification
PeptideIdentification.cpp # ConsensusMap
ProteinIdentification.cpp # MSExperiment / FileHandler

Acquisition.cpp
AcquisitionInfo.cpp
ContactPerson.cpp
DataArrays.cpp
Digestion.cpp
DocumentIDTagger.cpp
DocumentIdentifier.cpp

ChromatogramSettings.cpp
Precursor.cpp
PeptideEvidence.cpp

Identification.cpp
IdentificationHit.cpp
Instrument.cpp
InstrumentSettings.cpp

Modification.cpp
Tagging.cpp

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

