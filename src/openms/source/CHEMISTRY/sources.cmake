### the directory name
set(directory source/CHEMISTRY)

### list all filenames of the directory here
set(sources_list
ModificationDefinitionsSet.cpp # PeptideIdentification
NucleicAcidSpectrumGenerator.cpp # MSSpectrum
SpectrumAnnotator.cpp # MSSpectrum
SimpleTSGXLMS.cpp # OPXL
SvmTheoreticalSpectrumGenerator.cpp
SvmTheoreticalSpectrumGeneratorTrainer.cpp
SvmTheoreticalSpectrumGeneratorSet.cpp
TheoreticalSpectrumGenerator.cpp # MSSpectrum
TheoreticalSpectrumGeneratorXLMS.cpp # MSSpectrum
)

### add path to the filenames
set(sources)
foreach(i ${sources_list})
	list(APPEND sources ${directory}/${i})
endforeach(i)

### pass source file list to the upper instance
set(OpenMS_sources ${OpenMS_sources} ${sources})

### list all filenames of the directory here
set(sources_list
AASequence.cpp
CrossLinksDB.cpp
DecoyGenerator.cpp
DigestionEnzyme.cpp
DigestionEnzymeProtein.cpp
DigestionEnzymeRNA.cpp
DigestionEnzymeDB.cpp
EnzymaticDigestionLogModel.cpp
EnzymaticDigestion.cpp
ModificationDefinition.cpp
ModifiedNASequenceGenerator.cpp
ModifiedPeptideGenerator.cpp
Tagger.cpp


Element.cpp
ElementDB.cpp
EmpiricalFormula.cpp

ProteaseDB.cpp
ProteaseDigestion.cpp

NASequence.cpp
RNaseDB.cpp
RNaseDigestion.cpp
Ribonucleotide.cpp
RibonucleotideDB.cpp

Residue.cpp
ResidueDB.cpp
ResidueModification.cpp
ModificationsDB.cpp
WeightWrapper.cpp
)

### add path to the filenames
set(sources)
foreach(i ${sources_list})
	list(APPEND sources ${directory}/${i})
endforeach(i)

### pass source file list to the upper instance
set(OpenMS_BASE_sources ${OpenMS_BASE_sources} ${sources})

### source group definition
source_group("Source Files\\CHEMISTRY" FILES ${sources})
