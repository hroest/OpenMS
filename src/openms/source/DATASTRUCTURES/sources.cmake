### the directory name
set(directory source/DATASTRUCTURES)

set(OpenMS_CORE_sources ${OpenMS_CORE_sources} ${directory}/FASTAContainer.cpp)
set(OpenMS_MATH_sources ${OpenMS_MATH_sources} ${directory}/LPWrapper.cpp)

### list all filenames of the directory here
set(sources_list
DefaultParamHandler.cpp
QTCluster.cpp # AAsequence
GridFeature.cpp # BaseFeature/Peptide
ConstRefVector.cpp # FeatureMap
MatchedIterator.cpp # Math
OSWData.cpp # MSExperiment

Adduct.cpp # EmpiricalFormula
CalibrationData.cpp # Math
ChargePair.cpp # Adduct
Compomer.cpp # Adduct
MassExplainer.cpp # EmpiricalFormula

FlagSet.cpp
BinaryTreeNode.cpp
DataValue.cpp
Date.cpp
DateTime.cpp
Param.cpp
ParamValue.cpp
Map.cpp
Matrix.cpp
String.cpp
StringView.cpp
StringListUtils.cpp
StringUtils.cpp
StringUtilsSimple.cpp
StringConversions.cpp
ListUtils.cpp
ListUtilsIO.cpp
DIntervalBase.cpp
DPosition.cpp
DRange.cpp
DBoundingBox.cpp
ConvexHull2D.cpp
DistanceMatrix.cpp
CVMappingTerm.cpp
CVMappingRule.cpp
CVReference.cpp
CVMappings.cpp
ToolDescription.cpp
)

### add path to the filenames
set(sources)
foreach(i ${sources_list})
	list(APPEND sources ${directory}/${i})
endforeach(i)

### pass source file list to the upper instance
set(OpenMS_BASE_sources ${OpenMS_BASE_sources} ${sources})

### source group definition
source_group("Source Files\\DATASTRUCTURES" FILES ${sources})
