// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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

#include <OpenMS/KERNEL/OnDiscMSExperiment.h>

#include <OpenMS/FORMAT/MzMLFile.h>

namespace OpenMS
{

  OnDiscMSExperiment::OnDiscMSExperiment() 
  {
  }

  bool OnDiscMSExperiment::isSortedByRT() const
  {
    return meta_ms_experiment_->isSorted(false);
  }

  bool OnDiscMSExperiment::openFile(const String& filename, bool skipLoadingMetaData = false)
  {
    filename_ = filename;
    indexed_mzml_file_.openFile(filename);
    if (filename != "" && !skipLoadingMetaData)
    {
      loadMetaData_(filename);
    }
    return indexed_mzml_file_.getParsingSuccess();
  }

  void OnDiscMSExperiment::setSkipXMLChecks(bool skip)
  {
    indexed_mzml_file_.setSkipXMLChecks(skip);
  }

  OpenMS::Interfaces::SpectrumPtr OnDiscMSExperiment::getSpectrumById(Size id)
  {
    return indexed_mzml_file_.getSpectrumById(id);
  }

  OpenMS::Interfaces::ChromatogramPtr OnDiscMSExperiment::getChromatogramById(Size id)
  {
    return indexed_mzml_file_.getChromatogramById(id);
  }

  MSSpectrum OnDiscMSExperiment::getSpectrum(Size id)
  {
    if (id < meta_ms_experiment_->getNrSpectra())
    {
      MSSpectrum spectrum(meta_ms_experiment_->operator[](id));
      indexed_mzml_file_.getMSSpectrumById(static_cast<int>(id), spectrum);
      return spectrum;
    }
    else
    {
      MSSpectrum spectrum;
      indexed_mzml_file_.getMSSpectrumById(static_cast<int>(id), spectrum);
      return spectrum;
    }
  }

  MSChromatogram OnDiscMSExperiment::getChromatogram(Size id)
  {
    if (id < meta_ms_experiment_->getNrChromatograms())
    {
      MSChromatogram chromatogram(meta_ms_experiment_->getChromatogram(id));
      indexed_mzml_file_.getMSChromatogramById(static_cast<int>(id), chromatogram);
      return chromatogram;
    }
    else
    {
      MSChromatogram chromatogram;
      indexed_mzml_file_.getMSChromatogramById(static_cast<int>(id), chromatogram);
      return chromatogram;
    }
  }

  void OnDiscMSExperiment::loadMetaData_(const String& filename)
  {
    meta_ms_experiment_ = boost::shared_ptr< PeakMap >(new PeakMap);

    MzMLFile f;
    PeakFileOptions options = f.getOptions();
    options.setFillData(false);
    f.setOptions(options);
    f.load(filename, *meta_ms_experiment_.get());
  }

} //namespace OpenMS

