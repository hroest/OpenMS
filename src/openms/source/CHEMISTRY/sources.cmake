### the directory name
set(directory source/CHEMISTRY)

### list all filenames of the directory here
set(sources_list
AASequence.cpp
CrossLinksDB.cpp
DecoyGenerator.cpp
EnzymaticDigestionLogModel.cpp
EnzymaticDigestion.cpp
DigestionEnzyme.cpp
DigestionEnzymeProtein.cpp
DigestionEnzymeRNA.cpp
DigestionEnzymeDB.cpp
ModificationDefinition.cpp
ModificationDefinitionsSet.cpp
ModifiedNASequenceGenerator.cpp
ModifiedPeptideGenerator.cpp
NASequence.cpp
NucleicAcidSpectrumGenerator.cpp
ProteaseDB.cpp
ProteaseDigestion.cpp
RNaseDB.cpp
RNaseDigestion.cpp
Ribonucleotide.cpp
RibonucleotideDB.cpp
SpectrumAnnotator.cpp
SimpleTSGXLMS.cpp
SvmTheoreticalSpectrumGenerator.cpp
SvmTheoreticalSpectrumGeneratorTrainer.cpp
SvmTheoreticalSpectrumGeneratorSet.cpp
Tagger.cpp
TheoreticalSpectrumGenerator.cpp
TheoreticalSpectrumGeneratorXLMS.cpp
WeightWrapper.cpp
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
Residue.cpp
ResidueDB.cpp
ResidueModification.cpp
ModificationsDB.cpp

Element.cpp
ElementDB.cpp
EmpiricalFormula.cpp
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
