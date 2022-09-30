



    def getData(self):
        # cdef double[::1] arr = <double [:self.inst.get().data.size()]>self.inst.get().data.data()
        cdef shared_ptr[_OSBinaryDataArray] _r = self.inst
        cdef double[::1] arr = <double [:_r.get().data.size()]>_r.get().data.data()
        retval = np.asarray(arr)
        return retval
