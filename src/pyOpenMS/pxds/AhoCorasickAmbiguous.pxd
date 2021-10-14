## DELETE ?? 
from Types cimport *
from libcpp cimport bool
# from SeqanIncludeWrapper cimport *
from String cimport *

## TODO: do we need this here???
## kill it all!
## TODO: pulls in all of Seqan!

cdef extern from "<OpenMS/ANALYSIS/ID/AhoCorasickAmbiguous.h>" namespace "OpenMS":
    
    cdef cppclass AhoCorasickAmbiguous "OpenMS::AhoCorasickAmbiguous":
        AhoCorasickAmbiguous() nogil except +
        AhoCorasickAmbiguous(AhoCorasickAmbiguous) nogil except + #wrap-ignore
        # void initPattern(const PeptideDB & pep_db, const int aaa_max, const int mm_max, FuzzyACPattern & pattern) nogil except +
        AhoCorasickAmbiguous(const String & protein_sequence) nogil except +
        void setProtein(const String & protein_sequence) nogil except +
        # bool findNext(const FuzzyACPattern & pattern) nogil except +
        Size getHitDBIndex() nogil except +
        Int getHitProteinPosition() nogil except +

