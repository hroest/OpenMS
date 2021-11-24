### the directory name
set(directory source/DATASTRUCTURES)

### list all filenames of the directory here
set(sources_list
Adduct.cpp
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
MassExplainer.cpp
MatchedIterator.cpp
OSWData.cpp
QTCluster.cpp
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

### list all filenames of the directory here
set(sources_list
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
