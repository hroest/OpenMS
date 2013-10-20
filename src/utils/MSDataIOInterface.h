// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#ifndef OPENMS_LOWMEM_MSDATAIOINTERFACE
#define OPENMS_LOWMEM_MSDATAIOINTERFACE

#include <OpenMS/KERNEL/Peak1D.h>
#include <OpenMS/KERNEL/ChromatogramPeak.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/MSChromatogram.h>
#include <OpenMS/METADATA/ExperimentalSettings.h>

namespace OpenMS
{

  namespace Interfaces
  {

    /*
     * Provides the same interface as MSExperiment in order to be passed into mzMLHandler (for example).
     *
     * Only the interface used by the mzMLHandler is implemented here and this
     * class should be used as base class. The child class then implements those
     * functions which it actually expects to be used.
     *
     * Note that implementation of reset() is required.
    */
    template <typename PeakT = Peak1D, typename ChromatogramPeakT = ChromatogramPeak>
    class OPENMS_DLLAPI MSDataIOInterface :
      public ExperimentalSettings
    {

  public:

      /// @name Base type definitions
      //@{
      /// Peak type
      typedef PeakT PeakType;
      typedef MSSpectrum<PeakT> SpectrumType;
      /// Chromatogram type
      typedef ChromatogramPeakT ChromatogramPeakType;
      typedef MSChromatogram<ChromatogramPeakT> ChromatogramType;
      //@}

      virtual SpectrumType& operator[] (Size n)
      {
        throw Exception::NotImplemented(__FILE__, __LINE__, __PRETTY_FUNCTION__);
        return spectra_[n];
      }

      virtual const SpectrumType& operator[] (Size n) const
      {
        throw Exception::NotImplemented(__FILE__, __LINE__, __PRETTY_FUNCTION__);
        return spectra_[n];
      }

      virtual Size size() const
      {
        throw Exception::NotImplemented(__FILE__, __LINE__, __PRETTY_FUNCTION__);
        return spectra_.size(); 
      }

      virtual void reserve(Size /* s */)
      {
        throw Exception::NotImplemented(__FILE__, __LINE__, __PRETTY_FUNCTION__);
      }

      virtual bool empty() const
      {
        throw Exception::NotImplemented(__FILE__, __LINE__, __PRETTY_FUNCTION__);
        return spectra_.empty(); 
      }
      //@}

      virtual void reserveSpaceChromatograms(Size /* s */)
      {
        throw Exception::NotImplemented(__FILE__, __LINE__, __PRETTY_FUNCTION__);
      }

      virtual void reserveSpaceSpectra(Size /* s */)
      {
        throw Exception::NotImplemented(__FILE__, __LINE__, __PRETTY_FUNCTION__);
      }

      /// Constructor
      MSDataIOInterface() :
        ExperimentalSettings()
      {
      }

      /// Resets all internal values
      virtual void reset() = 0;

      /// adds a spectra to the list
      virtual void addSpectrum(MSSpectrum<PeakT> &)
      {
        throw Exception::NotImplemented(__FILE__, __LINE__, __PRETTY_FUNCTION__);
      }

      /// adds a chromatogram to the list
      virtual void addChromatogram(MSChromatogram<ChromatogramPeakType> &)
      {
        throw Exception::NotImplemented(__FILE__, __LINE__, __PRETTY_FUNCTION__);
      }

      virtual void resize(size_t) {}

      /// returns the spectra list
      virtual const std::vector<MSSpectrum<PeakT> > & getSpectra() const
      {
        throw Exception::NotImplemented(__FILE__, __LINE__, __PRETTY_FUNCTION__);
        return spectra_;
      }

      /// returns the spectra list
      virtual std::vector<MSSpectrum<PeakT> > & getSpectra() 
      {
        throw Exception::NotImplemented(__FILE__, __LINE__, __PRETTY_FUNCTION__);
        return spectra_;
      }

      /// sets the chromatogram list
      virtual void setChromatograms(const std::vector<MSChromatogram<ChromatogramPeakType> > & chromatograms)
      {
        throw Exception::NotImplemented(__FILE__, __LINE__, __PRETTY_FUNCTION__);
        chromatograms_ = chromatograms;
      }

      /// returns the chromatogram list
      virtual const std::vector<MSChromatogram<ChromatogramPeakType> > & getChromatograms() const
      {
        throw Exception::NotImplemented(__FILE__, __LINE__, __PRETTY_FUNCTION__);
        return chromatograms_;
      }

      /// returns a single chromatogram 
      virtual MSChromatogram<ChromatogramPeakType> & getChromatogram(Size id) 
      {
        throw Exception::NotImplemented(__FILE__, __LINE__, __PRETTY_FUNCTION__);
        return chromatograms_[id];
      }

      /**
        @brief Clears all data and meta data

        @param clear_meta_data If @em true, all meta data is cleared in addition to the data.
      */
      virtual void clear(bool)
      {
        throw Exception::NotImplemented(__FILE__, __LINE__, __PRETTY_FUNCTION__);
      }

  protected:

      /// chromatograms
      std::vector<MSChromatogram<ChromatogramPeakType> > chromatograms_;
      /// spectra
      std::vector<SpectrumType> spectra_;
    };

  } // Interfaces
} // NS OpenMS

#endif
