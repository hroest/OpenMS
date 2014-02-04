


cdef class ISpectrumAccess:

    cdef shared_ptr[_ISpectrumAccess] inst

    def __dealloc__(self):
         self.inst.reset()

    def init(self):
        assert False, 'Cannot initialize abstract base class'

    def getNrSpectra(self):
        assert self.inst.get() is not NULL, "There is no valid pointer for ISpectrumAccess"
        return self.inst.get().getNrSpectra()

    def getSpectrumById(self, id_):
        assert isinstance(id_, (int, long)), 'arg charge wrong type'
        assert self.inst.get() is not NULL, "There is no valid pointer for ISpectrumAccess"
    
        spec = OSSpectrum()
        spec.inst = self.inst.get().getSpectrumById((<int>id_))
        return spec 

    def getNrChromatograms(self):
        assert self.inst.get() is not NULL, "There is no valid pointer for ISpectrumAccess"
        return self.inst.get().getNrChromatograms()

    def getChromatogramById(self, id_):
        assert isinstance(id_, (int, long)), 'arg charge wrong type'
        assert self.inst.get() is not NULL, "There is no valid pointer for ISpectrumAccess"
    
        spec = OSChromatogram()
        spec.inst = self.inst.get().getChromatogramById((<int>id_))
        return spec 

