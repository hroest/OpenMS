from MSExperiment cimport *
from ChromatogramPeak cimport *
from Peak1D cimport *
from ProgressLogger cimport *
from TargetedExperiment cimport *
from TransformationDescription cimport *
from OpenSwathDataStructures cimport *
from ISpectrumAccess cimport *
from ChromatogramExtractorAlgorithm cimport *
from libcpp cimport bool

cdef extern from "<OpenMS/ANALYSIS/OPENSWATH/ChromatogramExtractor.h>" namespace "OpenMS":

    cdef cppclass ChromatogramExtractor(ProgressLogger):
        # wrap-inherits:
        #    ProgressLogger

        ChromatogramExtractor()                  nogil except +
        ChromatogramExtractor(ChromatogramExtractor)   nogil except + 

        void extractChromatograms(MSExperiment[Peak1D, ChromatogramPeak] & input,
                                  MSExperiment[Peak1D, ChromatogramPeak] & output, 
                                  TargetedExperiment & transition_exp,
                                  double extract_window,
                                  bool ppm,
                                  TransformationDescription trafo,
                                  double rt_extraction_window,
                                  String filter)

        ## wont work because of vector<shared_ptr<...> >
        ## void extractChromatograms(SpectrumAccessPtr input_, 
        ##     libcpp_vector[OSChromatogramPtr]& output, 
        ##     libcpp_vector[ExtractionCoordinates] extraction_coordinates,
        ##     double mz_extraction_window, bool ppm, String filter_) nogil except +

        # TODO immutable types by reference
        # void extract_value_tophat(MSSpectrum[Peak1D] input, double mz,
        #  Size peak_idx, double integrated_intensity, double extract_window, bool ppm)
        # void extract_value_bartlett(MSSpectrum[Peak1D] input, double mz, 
        #  Size peak_idx, double integrated_intensity, double extract_window, bool ppm)
    
