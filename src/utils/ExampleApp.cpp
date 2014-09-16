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
// $Maintainer: $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/config.h>

#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/DATASTRUCTURES/StringListUtils.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/KERNEL/ConversionHelper.h>

#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPickerHiRes.h>
#include <OpenMS/FILTERING/SMOOTHING/GaussFilter.h>

#include <OpenMS/FORMAT/DATAACCESS/MSDataWritingConsumer.h>
#include <OpenMS/FORMAT/DATAACCESS/MSDataChainingConsumer.h>
#include <OpenMS/INTERFACES/IMSDataConsumer.h>
#include <OpenMS/SYSTEM/SysInfo.h>
#include <OpenMS/KERNEL/OnDiscMSExperiment.h>
#include <OpenMS/FORMAT/IndexedMzMLFileLoader.h>

#include <OpenMS/INTERFACES/IMSDataConsumer.h>

using namespace OpenMS;
using namespace std;


/**
  @brief Helper class for the Low Memory peak-picking
*/
class PPHiResConsumer :
 public Interfaces::IMSDataConsumer< MSExperiment<> >
{

public:

  PPHiResConsumer(PeakPickerHiRes pp) :
      pp_(pp)
  {}

  void consumeSpectrum(MSSpectrum<>& s)
  {
    MSSpectrum<> sout;
    pp_.pick(s, sout);
    s = sout;
  }

  void consumeChromatogram(MSChromatogram<> & c) 
  {
    MSChromatogram<> cout;
    pp_.pick(c, cout);
    c = cout;
  }

  void setExpectedSize(Size expectedSpectra, Size expectedChromatograms) {}
  void setExperimentalSettings(const ExperimentalSettings& exp) {}

private:

  PeakPickerHiRes pp_;
};

class NoiseFilterGaussianConsumer :
 public Interfaces::IMSDataConsumer< MSExperiment<> >
{

public:

  NoiseFilterGaussianConsumer(GaussFilter nf) :
      nf_(nf)
  {}

  void consumeSpectrum(MSSpectrum<>& s)
  {
    nf_.filter(s);
  }

  void consumeChromatogram(MSChromatogram<> & c) 
  {
    nf_.filter(c);
  }

  void setExpectedSize(Size expectedSpectra, Size expectedChromatograms) {}
  void setExperimentalSettings(const ExperimentalSettings& exp) {}

private:

  GaussFilter nf_;
};

class ExampleApp :
  public TOPPBase
{
public:
  ExampleApp() :
    TOPPBase("ExampleApp", "")
  {
  }

protected:

  void registerOptionsAndFlags_()
  {
    registerInputFile_("in", "<file>", "", "Input file to convert.");
    setValidFormats_("in", ListUtils::create<String>("mzML"));

    registerOutputFile_("out", "<file>", "", "Output file");
    setValidFormats_("out", ListUtils::create<String>("mzML"));
    
    registerStringOption_("read_method", "<method>", "regular", "Method to read the file", false);
    String method("regular,indexed,streaming,cached");
    setValidStrings_("read_method", ListUtils::create<String>(method));
  }

  ExitCodes main_(int, const char**)
  {
    //-------------------------------------------------------------
    // parameter handling
    //-------------------------------------------------------------

    //input file names
    String in = getStringOption_("in");
    String out = getStringOption_("out");

    // Three steps: noise filter, peak picker, mzML writer
    GaussFilter nf;
    PeakPickerHiRes pp;
    NoiseFilterGaussianConsumer nf_consumer(nf);
    PPHiResConsumer pp_consumer(pp);
    PlainMSDataWritingConsumer wr_consumer(out);

    // Chain the steps
    MSDataChainingConsumer consumer;
    consumer.appendConsumer(&nf_consumer);
    consumer.appendConsumer(& pp_consumer);
    consumer.appendConsumer(& wr_consumer);

    // Transform file
    MzMLFile mzml;
    mzml.setLogType(log_type_);
    mzml.transform(in, &consumer, true);

    return EXECUTION_OK;
  }


};

int main(int argc, const char** argv)
{
  ExampleApp tool;
  return tool.main(argc, argv);
}

/// @endcond
