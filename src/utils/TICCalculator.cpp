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
// $Maintainer: Chris Bielow $
// $Authors: Marc Sturm, Andreas Bertsch, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/config.h>

#include <OpenMS/FORMAT/EDTAFile.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/FORMAT/IBSpectraFile.h>
#include <OpenMS/DATASTRUCTURES/StringListUtils.h>
#include <OpenMS/FORMAT/CachedMzML.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/KERNEL/ConversionHelper.h>

#include <OpenMS/FORMAT/DATAACCESS/MSDataWritingConsumer.h>
#include <OpenMS/SYSTEM/SysInfo.h>
#include <OpenMS/KERNEL/OnDiscMSExperiment.h>
#include <OpenMS/FORMAT/IndexedMzMLFileLoader.h>

#include <OpenMS/INTERFACES/IMSDataConsumer.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
  @page TOPP_FileConverter FileConverter

  @brief Converts between different MS file formats.

  <CENTER>
  <table>
  <tr>
  <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. predecessor tools </td>
  <td VALIGN="middle" ROWSPAN=3> \f$ \longrightarrow \f$ FileConverter \f$ \longrightarrow \f$</td>
  <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. successor tools </td>
  </tr>
  <tr>
  <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_GenericWrapper (e.g. for calling external converters) </td>
  <td VALIGN="middle" ALIGN = "center" ROWSPAN=2> any tool operating on the output format</td>
  </tr>
  <tr>
  <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> any vendor software exporting supported formats (e.g. mzXML) </td>
  </tr>
  </table>
  </CENTER>

  The main use of this tool is to convert data from external sources to the formats used by OpenMS/TOPP.
  Maybe most importantly, data from MS experiments in a number of different formats can be converted to mzML,
  the canonical file format used by OpenMS/TOPP for experimental data. (mzML is the PSI approved format and
  supports traceability of analysis steps.)

  Many different format conversions are supported, and some may be more useful than others. Depending on the
  file formats involved, information can be lost during conversion, e.g. when converting featureXML to mzData.
  In such cases a warning is shown.

  The input and output file types are determined from	the file extensions or from the first few lines of the
  files. If file type determination is not possible, the input or output file type has to be given explicitly.

  Conversion with the same output as input format is supported. In some cases, this can be helpful to remove
  errors from files, to update file formats to new versions, or to check whether information is lost upon
  reading or writing.

  Some information about the supported input types:
  @ref OpenMS::MzMLFile "mzML"
  @ref OpenMS::MzXMLFile "mzXML"
  @ref OpenMS::MzDataFile "mzData"
  @ref OpenMS::MascotGenericFile "mgf"
  @ref OpenMS::DTA2DFile "dta2d"
  @ref OpenMS::DTAFile "dta"
  @ref OpenMS::FeatureXMLFile "featureXML"
  @ref OpenMS::ConsensusXMLFile "consensusXML"
  @ref OpenMS::MS2File "ms2"
  @ref OpenMS::XMassFile "fid/XMASS"
  @ref OpenMS::MsInspectFile "tsv"
  @ref OpenMS::SpecArrayFile "peplist"
  @ref OpenMS::KroenikFile "kroenik"
  @ref OpenMS::EDTAFile "edta"

  See @ref TOPP_IDFileConverter for similar functionality for protein/peptide identification file formats.

  <B>The command line parameters of this tool are:</B>
  @verbinclude TOPP_FileConverter.cli
  <B>INI file documentation of this tool:</B>
  @htmlinclude TOPP_FileConverter.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES


class TICConsumer : 
    public Interfaces::IMSDataConsumer< MSExperiment<> >
{

    typedef MSExperiment<> MapType;
    typedef MapType::SpectrumType SpectrumType;
    typedef MapType::ChromatogramType ChromatogramType;

public:
  double TIC;
  int nr_spectra;
  long int nr_peaks;

  // Create new consumer, set TIC to zero
  TICConsumer() :
    TIC(0.0),
    nr_spectra(0.0)
    {}

  void consumeSpectrum(SpectrumType & s)
  {
    for (Size i = 0; i < s.size(); i++) 
    { 
      TIC += s[i].getIntensity(); 
    }
    nr_peaks += s.size();
    nr_spectra++;
  }

  void consumeChromatogram(ChromatogramType& /* c */) {}
  void setExpectedSize(Size expectedSpectra, Size expectedChromatograms) {}
  void setExperimentalSettings(const ExperimentalSettings& exp) {}
};

class TOPPFileConverter :
  public TOPPBase
{
public:
  TOPPFileConverter() :
    TOPPBase("FileConverter", "Converts between different MS file formats.")
  {
  }

protected:

  void registerOptionsAndFlags_()
  {
    registerInputFile_("in", "<file>", "", "Input file to convert.");
    registerStringOption_("in_type", "<type>", "", "Input file type -- default: determined from file extension or content\n", false);
    String formats("mzData,mzXML,mzML,cachedMzML,dta,dta2d,mgf,featureXML,consensusXML,ms2,fid,tsv,peplist,kroenik,edta");
    setValidFormats_("in", ListUtils::create<String>(formats));
    setValidStrings_("in_type", ListUtils::create<String>(formats));
    
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
    String read_method = getStringOption_("read_method");

    if (read_method == "streaming")
    {
      std::cout << "Read method: streaming" << std::endl;

      // Create the consumer, set output file name, transform
      TICConsumer consumer;
      MzMLFile mzml;
      mzml.setLogType(log_type_);
      mzml.transform(in, &consumer);

      std::cout << "There are " << consumer.nr_spectra << " spectra and " << consumer.nr_peaks << " peaks in the input file." << std::endl;
      std::cout << "The total ion current is " << consumer.TIC << std::endl;
      size_t after;
      SysInfo::getProcessMemoryConsumption(after);
      std::cout << " Memory consumption after " << after << std::endl;
    }
    else if (read_method == "regular")
    {
      std::cout << "Read method: regular" << std::endl;

      // convert MzXML to MzML
      MzMLFile mzml;
      mzml.setLogType(log_type_);
      MSExperiment<> map;
      mzml.load(in, map);
      double TIC = 0.0;
      long int nr_peaks = 0;
      for (Size i =0; i < map.size(); i++)
      {
        nr_peaks += map[i].size();
        for (Size j = 0; j < map[i].size(); j++)
        {
          TIC += map[i][j].getIntensity();
        }
      }

      std::cout << "There are " << map.size() << " spectra and " << nr_peaks << " peaks in the input file." << std::endl;
      std::cout << "The total ion current is " << TIC << std::endl;
      size_t after;
      SysInfo::getProcessMemoryConsumption(after);
      std::cout << " Memory consumption after " << after << std::endl;
    }
    else if (read_method == "indexed")
    {
      std::cout << "Read method: indexed" << std::endl;
      
      IndexedMzMLFileLoader imzml;
      //imzml.setLogType(log_type_);

      // load data from an indexed MzML file
      OnDiscMSExperiment<> map;
      imzml.load(in, map);
      // Get the first spectrum in memory, do some constant (non-changing) data processing
      double TIC = 0.0;
      long int nr_peaks = 0;
      for (Size i =0; i < map.getNrSpectra(); i++)
      {
        OpenMS::Interfaces::SpectrumPtr sptr = map.getSpectrumById(i);
        nr_peaks += sptr->getIntensityArray()->data.size();
        for (Size j = 0; j < sptr->getIntensityArray()->data.size(); j++)
        {
          TIC += sptr->getIntensityArray()->data[j];
        }
      }

      std::cout << "There are " << map.getNrSpectra() << " spectra and " << nr_peaks << " peaks in the input file." << std::endl;
      std::cout << "The total ion current is " << TIC << std::endl;
      size_t after;
      SysInfo::getProcessMemoryConsumption(after);
      std::cout << " Memory consumption after " << after << std::endl;
    }
    else if (read_method == "cached")
    {
      std::cout << "Read method: cached" << std::endl;


      // Special handling of cached mzML as input types: 
      // we expect two paired input files which we should read into exp
      std::vector<String> split_out;
      in.split(".cachedMzML", split_out);
      if (split_out.size() != 2)
      {
        LOG_ERROR << "Cannot deduce base path from input '" << in << "' (note that '.cachedMzML' should only occur once as the final ending)" << std::endl;
        return ILLEGAL_PARAMETERS;
      }
      String in_meta = split_out[0] + ".mzML";

      MzMLFile f;
      f.setLogType(log_type_);
      CachedmzML cacher;
      cacher.setLogType(log_type_);
      //MSExperiment<> tmp_exp;
      //MSExperiment<> exp;

      //f.load(in_meta, exp);
      // cacher.readMemdump(tmp_exp, in);

      CachedmzML cache;
      cache.createMemdumpIndex(in);
      const std::vector<std::streampos> spectra_index = cache.getSpectraIndex();
      const std::vector<std::streampos> chrom_index = cache.getChromatogramIndex();;

      std::ifstream ifs_;
      ifs_.open(in.c_str(), std::ios::binary);

      double TIC = 0.0;
      long int nr_peaks = 0;
      for (Size i=0; i < spectra_index.size(); ++i)
      {

        OpenSwath::BinaryDataArrayPtr mz_array(new OpenSwath::BinaryDataArray);
        OpenSwath::BinaryDataArrayPtr intensity_array(new OpenSwath::BinaryDataArray);
        int ms_level = -1;
        double rt = -1.0;
        ifs_.seekg(spectra_index[i]);
        CachedmzML::readSpectrumFast(mz_array, intensity_array, ifs_, ms_level, rt);

        nr_peaks += intensity_array->data.size();
        for (Size j = 0; j < intensity_array->data.size(); j++)
        {
          TIC += intensity_array->data[j];
        }
      }

      std::cout << "There are " << spectra_index.size() << " spectra and " << nr_peaks << " peaks in the input file." << std::endl;
      std::cout << "The total ion current is " << TIC << std::endl;
      size_t after;
      SysInfo::getProcessMemoryConsumption(after);
      std::cout << " Memory consumption after " << after << std::endl;
    }

  }


};

int main(int argc, const char** argv)
{
  TOPPFileConverter tool;
  size_t after, before;
  SysInfo::getProcessMemoryConsumption(before);
  std::cout << " Memory consumption before " << before << std::endl;
  tool.main(argc, argv);
  SysInfo::getProcessMemoryConsumption(after);
  std::cout << " Memory consumption after " << after << std::endl;
  // SysInfo::getProcessMemoryConsumption(before);
  return 0;
}

/// @endcond
