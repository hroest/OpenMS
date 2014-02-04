from Types cimport *
from libcpp.vector cimport vector as libcpp_vector
from OpenSwathDataStructures cimport *

cdef extern from "<OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/DATAACCESS/ISpectrumAccess.h>" namespace "OpenSwath":

  # This is an abstract base class in C++ with only pure virtual functions.
  # 
  # Thus we cannot create an instance of it but we can create a Python class holding a
  # pointer to a derived class of ISpectrumAccess.
  cdef cppclass ISpectrumAccess:
      # wrap-ignore

      shared_ptr[OSSpectrum] getSpectrumById(int id_) nogil except + # wrap-ignore
      size_t getNrSpectra() nogil except + # wrap-ignore

      # virtual std::vector<std::size_t> getSpectraByRT(double RT, double deltaRT) const = 0;
      # virtual SpectrumMeta getSpectrumMetaById(int id) const = 0;

      shared_ptr[OSChromatogram] getChromatogramById(int id_) nogil except + # wrap-ignore
      size_t getNrChromatograms() nogil except + # wrap-ignore
      libcpp_string getChromatogramNativeID(int id_) nogil except + # wrap-ignore

  ctypedef shared_ptr[ISpectrumAccess] SpectrumAccessPtr


