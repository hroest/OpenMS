from IsotopeDistribution cimport *

cdef extern from "<OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/IsotopePatternGenerator.h>" namespace "OpenMS":
    cdef cppclass IsotopePatternGenerator(IsotopeDistribution):

         # wrap-inherits:
         # IsotopeDistribution
         
         # wrap-ignore
         # ABSTRACT class
         # no-pxd-import

         IsotopePatternGenerator() # wrap-pass-constructor
         IsotopePatternGenerator(double probability_cutoff) nogil except +
         IsotopePatternGenerator(IsotopeDistribution) nogil except +
         # virtual void run(EmpiricalFormula&) = 0; # wrap-ignore
         void merge(double) nogil except +