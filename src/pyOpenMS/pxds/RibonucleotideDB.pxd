from Types cimport *
from String cimport *
from Ribonucleotide cimport *

cdef extern from "<OpenMS/CHEMISTRY/RibonucleotideDB.h>" namespace "OpenMS":
    
    cdef cppclass RibonucleotideDB "OpenMS::RibonucleotideDB":
        # wrap-manual-memory:
        #   cdef AutowrapPtrHolder[_RibonucleotideDB] inst

        RibonucleotideDB(RibonucleotideDB) nogil except + #wrap-ignore

        const Ribonucleotide * getRibonucleotide(String code) nogil except +
        const Ribonucleotide * getRibonucleotidePrefix(String code) nogil except +
        libcpp_pair[const Ribonucleotide *, const Ribonucleotide *] getRibonucleotideAlternatives(String code) nogil except + # wrap-ignore

# COMMENT: wrap static methods
cdef extern from "<OpenMS/CHEMISTRY/RibonucleotideDB.h>" namespace "OpenMS::RibonucleotideDB":
    
    RibonucleotideDB* getInstance() nogil except + # wrap-ignore

