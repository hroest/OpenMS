### the directory name
set(directory source/DATASTRUCTURES)

### list all filenames of the directory here
set(sources_list
Adduct.cpp
BinaryTreeNode.cpp
CalibrationData.cpp
ChargePair.cpp
Compomer.cpp
ConstRefVector.cpp
ConvexHull2D.cpp
CVMappingTerm.cpp
CVMappingRule.cpp
CVReference.cpp
CVMappings.cpp
DBoundingBox.cpp
DIntervalBase.cpp
DPosition.cpp
DRange.cpp
# DataValue.cpp
Date.cpp
DateTime.cpp
DefaultParamHandler.cpp
DistanceMatrix.cpp
FASTAContainer.cpp
FlagSet.cpp
GridFeature.cpp
#IsotopeCluster.h
#KDTree.h
ListUtils.cpp
ListUtilsIO.cpp
# LPWrapper.cpp
Map.cpp
MassExplainer.cpp
MatchedIterator.cpp
Matrix.cpp
OSWData.cpp
Param.cpp
# ParamValue.cpp
QTCluster.cpp
# String.cpp
StringView.cpp
StringListUtils.cpp
StringUtils.cpp
StringUtilsSimple.cpp
StringConversions.cpp
ToolDescription.cpp
)

### add path to the filenames
set(sources)
foreach(i ${sources_list})
	list(APPEND sources ${directory}/${i})
endforeach(i)

### pass source file list to the upper instance
set(OpenMS_CORE_sources ${OpenMS_CORE_sources} ${sources})

set(OpenMS_MATH_sources ${OpenMS_MATH_sources} ${directory}/LPWrapper.cpp)

set(OpenMS_BASE_sources ${OpenMS_BASE_sources} ${directory}/String.cpp)
set(OpenMS_BASE_sources ${OpenMS_BASE_sources} ${directory}/DataValue.cpp)
set(OpenMS_BASE_sources ${OpenMS_BASE_sources} ${directory}/ParamValue.cpp)

### source group definition
source_group("Source Files\\DATASTRUCTURES" FILES ${sources})
