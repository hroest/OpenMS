from Types cimport *
from libcpp cimport bool
from Types cimport *
from String cimport *
from InterfaceDataStructures cimport *
# from ISpectrumAccess cimport *
from IndexedMzMLDecoder cimport *
from MSSpectrum cimport *
from MSChromatogram cimport *

cdef extern from "<OpenMS/FORMAT/HANDLERS/IndexedMzMLHandler.h>" namespace "OpenMS":
    
    cdef cppclass Internal_IndexedMzMLHandler "OpenMS::Internal::IndexedMzMLHandler":

        Internal_IndexedMzMLHandler() nogil except +
        Internal_IndexedMzMLHandler(Internal_IndexedMzMLHandler) nogil except +
        Internal_IndexedMzMLHandler(String filename) nogil except +

        void openFile(String filename) nogil except +
        bool getParsingSuccess() nogil except +
        size_t getNrSpectra() nogil except +
        size_t getNrChromatograms() nogil except +
        shared_ptr[Spectrum] getSpectrumById(int id_) nogil except +
        shared_ptr[Chromatogram] getChromatogramById(int id_) nogil except +
        MSSpectrum getMSSpectrumById(int id_) nogil except +
        MSChromatogram getMSChromatogramById(int id_) nogil except +
        void setSkipXMLChecks(bool skip) nogil except +

