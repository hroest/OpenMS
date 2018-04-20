from Types cimport *

from PeakFileOptions cimport *
from MzMLFile cimport *
from OnDiscMSExperiment cimport *

cdef extern from "<OpenMS/FORMAT/IndexedMzMLFile.h>" namespace "OpenMS":

    cdef cppclass IndexedMzMLFile:

        IndexedMzMLFile() nogil except +
 
        bool load(String, OnDiscMSExperiment &) nogil except+
        void store(String, OnDiscMSExperiment &) nogil except+
        void store(String, MSExperiment &) nogil except+

        PeakFileOptions getOptions() nogil except +
        void setOptions(PeakFileOptions) nogil except +

