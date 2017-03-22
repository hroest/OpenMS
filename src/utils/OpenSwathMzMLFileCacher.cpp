// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2016.
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

#include "OpenMS/FORMAT/CachedMzML.h"

#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/FORMAT/MzMLFile.h>

#include <fstream>


using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
  @page UTILS_OpenSwathMzMLFileCacher OpenSwathMzMLFileCacher

  @brief Serialize a spectra and/or chromatogram mzML file

  This class will serialize a spectra and/or chromatogram mzML file and store
  it in a binary format that contains ONLY the spectra and chromatogram data
  (no metadata).
 
  This is implemented using the write_memdump and read_memdump functions.
  For reading there are 2 options
  - read the whole file into the OpenMS datastructures
  - read only an index (read_memdump_idx) of the spectra and chromatograms and then use
    random-access to retrieve a specific spectra from the disk (read_memdump_spectra)

  @note This tool is experimental!

  <B>The command line parameters of this tool are:</B>
  @verbinclude UTILS_OpenSwathMzMLFileCacher.cli
  <B>INI file documentation of this tool:</B>
  @htmlinclude UTILS_OpenSwathMzMLFileCacher.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

#include <sqlite3.h>
#include <zlib.h>
#include <OpenMS/FORMAT/Base64.h>
#include <OpenMS/FORMAT/MSNumpressCoder.h>
#include <QByteArray>
class OPENMS_DLLAPI MSNumpressCoder_Internal : MSNumpressCoder
{

  public :

    MSNumpressCoder_Internal() {}

    void encodeNP_raw(const std::vector<double> & in, String & result, const MSNumpressCoder::NumpressConfig & config)
    {
      MSNumpressCoder::encodeNP_(in, result, config);
    }

    void decodeNP_raw(const std::string & in, std::vector<double> & out,
        const NumpressConfig & config)
    {
      decodeNP_(in, out, config);
    }

};

namespace OpenMS
{
  namespace Internal
  {

    static void compress_str(std::string& str, std::string& compressed)
    {
      compressed.clear();

      unsigned long sourceLen =   (unsigned long)str.size();
      unsigned long compressed_length = //compressBound((unsigned long)str.size());
                                        sourceLen + (sourceLen >> 12) + (sourceLen >> 14) + 11; // taken from zlib's compress.c, as we cannot use compressBound*

      int zlib_error;
      do
      {
        compressed.resize(compressed_length);
        zlib_error = compress(reinterpret_cast<Bytef*>(&compressed[0]), &compressed_length, reinterpret_cast<Bytef*>(&str[0]), (unsigned long) str.size());

        switch (zlib_error)
        {
        case Z_MEM_ERROR:
          throw Exception::OutOfMemory(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, compressed_length);

        case Z_BUF_ERROR:
          compressed_length *= 2;
        }
      }
      while (zlib_error == Z_BUF_ERROR);

      if (zlib_error != Z_OK)
      {
        throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Compression error?");
      }
      compressed.resize(compressed_length);
    }

    static void uncompress_str(const void * tt, size_t blob_bytes, std::string& uncompressed)
    {
      // take a leap of faith and assume the input is valid
      uncompressed.clear();
      QByteArray raw_data = QByteArray::fromRawData((const char*)tt, blob_bytes);

      // base64_uncompressed = QByteArray::fromBase64(herewego);
      // void Base64::decodeSingleString(const String& in, QByteArray& base64_uncompressed, bool zlib_compression)
      // base64_uncompressed = QByteArray::fromBase64(herewego);
      // QByteArray base64_uncompressed = str.c_str();
      QByteArray base64_uncompressed = raw_data;
      if (true)
      {
        QByteArray czip;
        czip.resize(4);
        czip[0] = (base64_uncompressed.size() & 0xff000000) >> 24;
        czip[1] = (base64_uncompressed.size() & 0x00ff0000) >> 16;
        czip[2] = (base64_uncompressed.size() & 0x0000ff00) >> 8;
        czip[3] = (base64_uncompressed.size() & 0x000000ff);
        czip += base64_uncompressed;
        base64_uncompressed = qUncompress(czip);

        if (base64_uncompressed.isEmpty())
        {
          throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Decompression error?");
        }
      }

      // Note that we may have zero bytes in the string, so we cannot use QString
      uncompressed = std::string(base64_uncompressed.data(), base64_uncompressed.size());
    }

  }
}

class OPENMS_DLLAPI SqMassWriter
{
  String filename_;

  /// Decoder/Encoder for Base64-data in MzML
  Base64 base64coder_;
  MSNumpressCoder numpress_coder_;

  Int spec_id_;
  Int chrom_id_;

  public:

    SqMassWriter(String filename) :
      filename_(filename),
      base64coder_(),
      numpress_coder_(),
      spec_id_(0),
      chrom_id_(0)
    {
    }

    static int callback(void *NotUsed, int argc, char **argv, char **azColName)
    {
      int i;
      for (i=0; i<argc; i++)
      {
        printf("%s = %s\n", azColName[i], argv[i] ? argv[i] : "NULL");
      }
      printf("\n");
      return(0);
    }

    void init()
    {
      sqlite3 *db;
      char *zErrMsg = 0;
      int rc;
      char *create_sql;

      // delete file if present
      // TODO
      // remove(filename_);

      // Open database
      rc = sqlite3_open(filename_.c_str(), &db);
      if (rc)
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__, String("Can't open database: ") + sqlite3_errmsg(db));
      }

      // Create SQL structure
      create_sql = 

    // data table
    //  - compression is one of 0 = no, 1 = zlib, 2 = np-linear, 3 = np-slof, 4 = np-pic, 5 = np-linear + zlib, 6 = np-slof + zlib, 7 = np-pic + zlib
    //  - data_type is one of 0 = mz, 1 = int, 2 = rt
    //  - data contains the raw (blob) data for a single data array
            "CREATE TABLE DATA(" \
            "SPECTRUM_ID INT," \
            "CHROMATOGRAM_ID INT," \
            "COMPRESSION INT," \
            "DATA_TYPE INT," \
            "DATA BLOB NOT NULL" \
            ");" \

    // spectrum table
            "CREATE TABLE SPECTRUM(" \
            "ID INT PRIMARY KEY NOT NULL," \
            "MSLEVEL INT NULL," \
            "RETENTION_TIME REAL NULL," \
            "NATIVE_ID TEXT NOT NULL" \
            ");" \

    // chromatogram table
            "CREATE TABLE CHROMATOGRAM(" \
            "ID INT PRIMARY KEY NOT NULL," \
            "NATIVE_ID TEXT NOT NULL" \
            ");" \
    
    // product table
            "CREATE TABLE PRODUCT(" \
            "SPECTRUM_ID INT," \
            "CHROMATOGRAM_ID INT," \
            "CHARGE INT NULL," \
            "ISOLATION_TARGET REAL NULL," \
            "ISOLATION_LOWER REAL NULL," \
            "ISOLATION_UPPER REAL NULL" \
            ");" \

    // precursor table
            "CREATE TABLE PRECURSOR(" \
            "SPECTRUM_ID INT," \
            "CHROMATOGRAM_ID INT," \
            "CHARGE INT NULL," \
            "PEPTIDE_SEQUENCE TEXT NULL," \
            "DRIFT_TIME REAL NULL," \
            "ISOLATION_TARGET REAL NULL," \
            "ISOLATION_LOWER REAL NULL," \
            "ISOLATION_UPPER REAL NULL" \
            ");";


      // Execute SQL statement
      rc = sqlite3_exec(db, create_sql, callback, 0, &zErrMsg);
      if( rc != SQLITE_OK )
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__, 
            zErrMsg);
        sqlite3_free(zErrMsg);
      }
      else {
        std::cout << "Done creating tables" << std::endl;
      }
      sqlite3_close(db);
    }

    void executeBlobBind(sqlite3 *db, String& prepare_statement, std::vector<String>& data)
    {
      int rc;

      // The calling procedure is responsible for deleting the compiled SQL statement using sqlite3_finalize() after it has finished with it.
      sqlite3_stmt *stmt = NULL;
      const char *curr_loc;
      rc = sqlite3_prepare_v2(db, prepare_statement.c_str(), prepare_statement.size(), &stmt, &curr_loc);
      if (rc != SQLITE_OK) { throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__, sqlite3_errmsg(db)); }

      for (Size k = 0; k < data.size(); k++)
      {
        // Fifth argument is a destructor for the blob.
        // SQLITE_STATIC because the statement is finalized
        // before the buffer is freed:
        rc = sqlite3_bind_blob(stmt, k+1, data[k].c_str(), data[k].size(), SQLITE_STATIC);
        if (rc != SQLITE_OK) { throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__, sqlite3_errmsg(db)); } 
      }

      rc = sqlite3_step(stmt);
      if (rc != SQLITE_DONE) { throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__, sqlite3_errmsg(db)); }

      // free memory again
      sqlite3_finalize(stmt);
    }

    void executeSql(sqlite3 *db, const std::stringstream& statement)
    {
      char *zErrMsg = 0;
      std::string insert_str = statement.str();
      int rc = sqlite3_exec(db, insert_str.c_str(), callback, 0, &zErrMsg);
      if( rc != SQLITE_OK )
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__, 
            zErrMsg);
        sqlite3_free(zErrMsg);
      }
    }

    void writeSpectra(const std::vector<MSSpectrum<> >& spectra)
    {
      if (spectra.empty()) return;

      sqlite3 *db;
      char *zErrMsg = 0;
      int rc;

      // Open database
      rc = sqlite3_open(filename_.c_str(), &db);
      if (rc)
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__, String("Can't open database: ") + sqlite3_errmsg(db));
      }

      // prepare streams and set required precision (default is 6 digits)
      std::stringstream insert_spectra_sql;
      std::stringstream insert_precursor_sql;
      std::stringstream insert_product_sql;

      insert_spectra_sql.precision(11);
      insert_precursor_sql.precision(11);
      insert_product_sql.precision(11);

      // Encoding options
      MSNumpressCoder::NumpressConfig npconfig_mz;
      npconfig_mz.estimate_fixed_point = true; // critical
      npconfig_mz.numpressErrorTolerance = -1.0; // skip check, faster
      npconfig_mz.setCompression("linear");
      npconfig_mz.linear_fp_mass_acc = 0.0001; // set the desired mass accuracy = 1ppm at 100 m/z
      MSNumpressCoder::NumpressConfig npconfig_int;
      npconfig_int.estimate_fixed_point = true; // critical
      npconfig_int.numpressErrorTolerance = -1.0; // skip check, faster
      npconfig_int.setCompression("slof");

      String prepare_statement = "INSERT INTO DATA(SPECTRUM_ID, DATA_TYPE, COMPRESSION, DATA) VALUES ";
      std::vector<String> data;
      int sql_it = 1;
      int nr_precursors = 0;
      int nr_products = 0;
      for (Size k = 0; k < spectra.size(); k++)
      {
        const MSSpectrum<>& spec = spectra[k];
        insert_spectra_sql << "INSERT INTO SPECTRUM(ID, NATIVE_ID, MSLEVEL, RETENTION_TIME) VALUES (" << spec_id_ << ",'" << spec.getNativeID() << "'," 
          << spec.getMSLevel() << "," << spec.getRT() << "); ";

        if (!spec.getPrecursors().empty())
        {
          if (spec.getPrecursors().size() > 1) std::cout << "WARNING cannot store more than first precursor" << std::endl;
          OpenMS::Precursor prec = spec.getPrecursors()[0];
          String pepseq;
          if (prec.metaValueExists("peptide_sequence"))
          {
            pepseq = prec.getMetaValue("peptide_sequence");
            insert_precursor_sql << "INSERT INTO PRECURSOR (SPECTRUM_ID, CHARGE, ISOLATION_TARGET, PEPTIDE_SEQUENCE) VALUES (" << spec_id_ << "," << 
              prec.getCharge() << "," << prec.getMZ() << ",'" << pepseq << "'" <<  "); ";
          }
          else
          {
            insert_precursor_sql << "INSERT INTO PRECURSOR (SPECTRUM_ID, CHARGE, ISOLATION_TARGET) VALUES (" << spec_id_ << "," << prec.getCharge() << "," << prec.getMZ() << "); ";
          }
          nr_products++;
        }

        if (!spec.getProducts().empty())
        {
          if (spec.getProducts().size() > 1) std::cout << "WARNING cannot store more than first product" << std::endl;
          OpenMS::Product prod = spec.getProducts()[0];
          insert_product_sql << "INSERT INTO PRODUCT (SPECTRUM_ID, CHARGE, ISOLATION_TARGET) VALUES (" << spec_id_ << "," << 0 << "," << prod.getMZ() <<  "); ";
          nr_products++;
        }

        //  data_type is one of 0 = mz, 1 = int, 2 = rt
        //  compression is one of 0 = no, 1 = zlib, 2 = np-linear, 3 = np-slof, 4 = np-pic, 5 = np-linear + zlib, 6 = np-slof + zlib, 7 = np-pic + zlib

        // encode mz data
        {
          std::vector<double> data_to_encode;
          data_to_encode.resize(spec.size());
          for (Size p = 0; p < spec.size(); ++p)
          {
            data_to_encode[p] = spec[p].getMZ();
          }

          String uncompressed_str;
          String encoded_string;
          MSNumpressCoder_Internal().encodeNP_raw(data_to_encode, uncompressed_str, npconfig_mz);
          OpenMS::Internal::compress_str(uncompressed_str, encoded_string);
          data.push_back(encoded_string);
          prepare_statement += String("(") + spec_id_ + ", 0, 5, ?" + sql_it++ + " ),";
        }

        // encode intensity data
        {
          std::vector<double> data_to_encode;
          data_to_encode.resize(spec.size());
          for (Size p = 0; p < spec.size(); ++p)
          {
            data_to_encode[p] = spec[p].getIntensity();
          }

          String uncompressed_str;
          String encoded_string;
          MSNumpressCoder_Internal().encodeNP_raw(data_to_encode, uncompressed_str, npconfig_int);
          OpenMS::Internal::compress_str(uncompressed_str, encoded_string);
          data.push_back( encoded_string );
          prepare_statement += String("(") + spec_id_ + ", 1, 6, ?" + sql_it++ + " ),";
        }
        spec_id_++;

        if (sql_it > 500) // flush after 500 spectra as sqlite can only handle so many bind_blob statments
        {
          prepare_statement.resize( prepare_statement.size() -1 ); // remove last ","
          executeBlobBind(db, prepare_statement, data);

          data.clear();
          prepare_statement = "INSERT INTO DATA(SPECTRUm_ID, DATA_TYPE, COMPRESSION, DATA) VALUES ";
          sql_it = 1;
        }

      }

      prepare_statement.resize( prepare_statement.size() -1 );
      executeBlobBind(db, prepare_statement, data);

      sqlite3_exec(db, "BEGIN TRANSACTION", NULL, NULL, &zErrMsg);

      executeSql(db, insert_spectra_sql);
      if (nr_precursors > 0) executeSql(db, insert_precursor_sql);
      if (nr_products > 0) executeSql(db, insert_product_sql);

      sqlite3_exec(db, "END TRANSACTION", NULL, NULL, &zErrMsg);

      sqlite3_close(db);
    }

    void writeChromatograms(const std::vector<MSChromatogram<> >& chroms)
    {
      if (chroms.empty()) return;

      sqlite3 *db;
      char *zErrMsg = 0;
      int rc;

      // Open database
      rc = sqlite3_open(filename_.c_str(), &db);
      if (rc)
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__, String("Can't open database: ") + sqlite3_errmsg(db));
      }

      // prepare streams and set required precision (default is 6 digits)
      std::stringstream insert_chrom_sql;
      std::stringstream insert_precursor_sql;
      std::stringstream insert_product_sql;

      insert_chrom_sql.precision(11);
      insert_precursor_sql.precision(11);
      insert_product_sql.precision(11);

      // Encoding options
      MSNumpressCoder::NumpressConfig npconfig_mz;
      npconfig_mz.estimate_fixed_point = true; // critical
      npconfig_mz.numpressErrorTolerance = -1.0; // skip check, faster
      npconfig_mz.setCompression("linear");
      npconfig_mz.linear_fp_mass_acc = 0.05; // set the desired RT accuracy (0.05 seconds)
      MSNumpressCoder::NumpressConfig npconfig_int;
      npconfig_int.estimate_fixed_point = true; // critical
      npconfig_int.numpressErrorTolerance = -1.0; // skip check, faster
      npconfig_int.setCompression("slof");

      String prepare_statement = "INSERT INTO DATA(CHROMATOGRAM_ID, DATA_TYPE, COMPRESSION, DATA) VALUES ";
      std::vector<String> data;
      int sql_it = 1;
      for (Size k = 0; k < chroms.size(); k++)
      {
        const MSChromatogram<>& chrom = chroms[k];
        insert_chrom_sql << "INSERT INTO CHROMATOGRAM (ID, NATIVE_ID) VALUES (" << chrom_id_ << ",'" << chrom.getNativeID() << "'); ";

        OpenMS::Precursor prec = chrom.getPrecursor();
        String pepseq;
        if (prec.metaValueExists("peptide_sequence"))
        {
          pepseq = prec.getMetaValue("peptide_sequence");
          insert_precursor_sql << "INSERT INTO PRECURSOR (CHROMATOGRAM_ID, CHARGE, ISOLATION_TARGET, PEPTIDE_SEQUENCE) VALUES (" << chrom_id_ << "," << 
            prec.getCharge() << "," << prec.getMZ() << ",'" << pepseq << "'" <<  "); ";
        }
        else
        {
          insert_precursor_sql << "INSERT INTO PRECURSOR (CHROMATOGRAM_ID, CHARGE, ISOLATION_TARGET) VALUES (" << chrom_id_ << "," << prec.getCharge() << "," << prec.getMZ() << "); ";
        }

        OpenMS::Product prod = chrom.getProduct();
        insert_product_sql << "INSERT INTO PRODUCT (CHROMATOGRAM_ID, CHARGE, ISOLATION_TARGET) VALUES (" << chrom_id_ << "," << 0 << "," << prod.getMZ() <<  "); ";

        //  data_type is one of 0 = mz, 1 = int, 2 = rt
        //  compression is one of 0 = no, 1 = zlib, 2 = np-linear, 3 = np-slof, 4 = np-pic, 5 = np-linear + zlib, 6 = np-slof + zlib, 7 = np-pic + zlib

        // encode retention time data
        {
          std::vector<double> data_to_encode;
          data_to_encode.resize(chrom.size());
          for (Size p = 0; p < chrom.size(); ++p)
          {
            data_to_encode[p] = chrom[p].getRT();
          }

          String uncompressed_str;
          String encoded_string;
          MSNumpressCoder_Internal().encodeNP_raw(data_to_encode, uncompressed_str, npconfig_mz);
          OpenMS::Internal::compress_str(uncompressed_str, encoded_string);
          data.push_back(encoded_string);
          prepare_statement += String("(") + chrom_id_ + ", 2, 5, ?" + sql_it++ + " ),";
        }

        // encode intensity data
        {
          std::vector<double> data_to_encode;
          data_to_encode.resize(chrom.size());
          for (Size p = 0; p < chrom.size(); ++p)
          {
            data_to_encode[p] = chrom[p].getIntensity();
          }

          String uncompressed_str;
          String encoded_string;
          MSNumpressCoder_Internal().encodeNP_raw(data_to_encode, uncompressed_str, npconfig_int);
          OpenMS::Internal::compress_str(uncompressed_str, encoded_string);
          data.push_back( encoded_string );
          prepare_statement += String("(") + chrom_id_ + ", 1, 6, ?" + sql_it++ + " ),";
        }
        chrom_id_++;

        if (sql_it > 500) // flush after 500 chromatograms as sqlite can only handle so many bind_blob statments
        {
          prepare_statement.resize( prepare_statement.size() -1 ); // remove last ","
          executeBlobBind(db, prepare_statement, data);

          data.clear();
          prepare_statement = "INSERT INTO DATA(CHROMATOGRAM_ID, DATA_TYPE, COMPRESSION, DATA) VALUES ";
          sql_it = 1;
        }

      }

      prepare_statement.resize( prepare_statement.size() -1 );
      executeBlobBind(db, prepare_statement, data);

      sqlite3_exec(db, "BEGIN TRANSACTION", NULL, NULL, &zErrMsg);

      executeSql(db, insert_chrom_sql);
      executeSql(db, insert_precursor_sql);
      executeSql(db, insert_product_sql);

      sqlite3_exec(db, "END TRANSACTION", NULL, NULL, &zErrMsg);

      sqlite3_close(db);
    }

};

class OPENMS_DLLAPI SqMassReader
{

  String filename_; 

  /// Decoder/Encoder for Base64-data in MzML
  Base64 base64coder_;
  MSNumpressCoder numpress_coder_;

  public:

    SqMassReader(String filename) :
      filename_(filename),
      base64coder_(),
      numpress_coder_()
    {
    }

    static int callback(void * /* NotUsed */, int argc, char **argv, char **azColName)
    {
      int i;
      for (i=0; i<argc; i++)
      {
        printf("%s = %s\n", azColName[i], argv[i] ? argv[i] : "NULL");
      }
      printf("\n");
      return(0);
    }

    void read(MSExperiment<>& exp, bool meta_only = false)
    {
      sqlite3 *db;
      sqlite3_stmt * stmt;
      // char *zErrMsg = 0;
      int rc;
      std::string select_sql;
      // const char* data = "Callback function called";

      // Open database
      rc = sqlite3_open(filename_.c_str(), &db);
      if (rc)
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__, String("Can't open database: ") + sqlite3_errmsg(db));
      }

      // creates the chromatograms but does not fill them with data (provides option to return meta-data only)
      std::vector<MSChromatogram<> > chromatograms;
      prepareChroms_(db, chromatograms);
      if (meta_only) 
      {
        exp.setChromatograms(chromatograms);
        return;
      }

      select_sql = "SELECT " \
                   "CHROMATOGRAM.ID as chrom_id," \
                   "CHROMATOGRAM.NATIVE_ID as chrom_native_id," \
                   "DATA.COMPRESSION as data_compression," \
                   "DATA.DATA_TYPE as data_type," \
                   "DATA.DATA as binary_data " \
                   "FROM CHROMATOGRAM " \
                   "INNER JOIN DATA ON CHROMATOGRAM.ID = DATA.CHROMATOGRAM_ID " \
                   ";";


      /* Execute SQL statement */
      sqlite3_prepare(db, select_sql.c_str(), -1, &stmt, NULL);
      sqlite3_step( stmt );

      // TODO ensure that all chroms have their data...
      std::vector<int> chromdata; chromdata.resize(chromatograms.size());
      int k = 0;
      while (sqlite3_column_type( stmt, 0 ) != SQLITE_NULL)
      {
        int chrom_id = sqlite3_column_int( stmt, 0 );
        const unsigned char * native_id_ = sqlite3_column_text(stmt, 1);
        std::string native_id(reinterpret_cast<const char*>(native_id_), sqlite3_column_bytes(stmt, 1));

        if (chrom_id >= chromatograms.size())
        {
          throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__, 
              "Data for non-existent chromatogram found");
        }
        if (native_id != chromatograms[chrom_id].getNativeID())
        {
          throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__, 
              "Native id for chromatogram doesnt match");
        }

        int compression = sqlite3_column_int( stmt, 2 );
        int data_type = sqlite3_column_int( stmt, 3 );

        const void * tt = sqlite3_column_blob(stmt, 4);
        size_t blob_bytes = sqlite3_column_bytes(stmt, 4);

        // data_type is one of 0 = mz, 1 = int, 2 = rt
        // compression is one of 0 = no, 1 = zlib, 2 = np-linear, 3 = np-slof, 4 = np-pic, 5 = np-linear + zlib, 6 = np-slof + zlib, 7 = np-pic + zlib
        std::vector<double> data;
        if (compression == 5)
        {
          std::string uncompressed;
          OpenMS::Internal::uncompress_str(tt, blob_bytes, uncompressed);
          MSNumpressCoder::NumpressConfig config;
          config.setCompression("linear");
          MSNumpressCoder_Internal().decodeNP_raw(uncompressed, data, config);
        }
        else if (compression == 6)
        {
          std::string uncompressed;
          OpenMS::Internal::uncompress_str(tt, blob_bytes, uncompressed);
          MSNumpressCoder::NumpressConfig config;
          config.setCompression("slof");
          MSNumpressCoder_Internal().decodeNP_raw(uncompressed, data, config);
        }
        else
        {
          throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__, 
              "Compression not supported");
        }

        if (data_type == 1)
        {
          // intensity
          if (chromatograms[chrom_id].empty()) chromatograms[chrom_id].resize(data.size());
          std::vector< double >::iterator data_it = data.begin();
          for (MSChromatogram<>::iterator it = chromatograms[chrom_id].begin(); it != chromatograms[chrom_id].end(); it++, data_it++)
          {
            it->setIntensity(*data_it);
          }
          chromdata[chrom_id] += 1;
        }
        else if (data_type == 2)
        {
          // rt
          if (chromatograms[chrom_id].empty()) chromatograms[chrom_id].resize(data.size());
          std::vector< double >::iterator data_it = data.begin();
          for (MSChromatogram<>::iterator it = chromatograms[chrom_id].begin(); it != chromatograms[chrom_id].end(); it++, data_it++)
          {
            it->setRT(*data_it);
          }
          chromdata[chrom_id] += 1;
        }
        else
        {
          throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__, 
              "Found data type other than RT/Intensity for chromatograms");
        }

        sqlite3_step( stmt );
        k++;
      }

      sqlite3_finalize(stmt);
    
      exp.setChromatograms(chromatograms);
    }

  private:

    void prepareChroms_(sqlite3 *db, std::vector<MSChromatogram<> >& chromatograms)
    {
      sqlite3_stmt * stmt;
      std::string select_sql;
      select_sql = "SELECT " \
                   "CHROMATOGRAM.ID as chrom_id," \
                   "CHROMATOGRAM.NATIVE_ID as chrom_native_id," \
                   "PRECURSOR.CHARGE as precursor_charge," \
                   "PRECURSOR.DRIFT_TIME as precursor_dt," \
                   "PRECURSOR.ISOLATION_TARGET as precursor_mz," \
                   "PRECURSOR.ISOLATION_LOWER as precursor_mz_lower," \
                   "PRECURSOR.ISOLATION_UPPER as precursor_mz_upper," \
                   "PRECURSOR.PEPTIDE_SEQUENCE as precursor_seq," \
                   "PRODUCT.CHARGE as product_charge," \
                   "PRODUCT.ISOLATION_TARGET as product_mz," \
                   "PRODUCT.ISOLATION_LOWER as product_mz_lower," \
                   "PRODUCT.ISOLATION_UPPER as product_mz_upper " \
                   "FROM CHROMATOGRAM " \
                   "INNER JOIN PRECURSOR ON CHROMATOGRAM.ID = PRECURSOR.CHROMATOGRAM_ID " \
                   "INNER JOIN PRODUCT ON CHROMATOGRAM.ID = PRODUCT.CHROMATOGRAM_ID " \
                   ";";

      /// TODO : do we want to support reading a subset of the data (e.g. only chromatograms xx - yy)
      ///   readChromatograms_(db, stmt, chromatograms);
      /// }

      /// void readChromatograms_(sqlite3 *db, sqlite3_stmt* stmt, std::vector<MSChromatogram<> >& chromatograms)
      /// {

      // See https://www.sqlite.org/c3ref/column_blob.html
      // The pointers returned are valid until a type conversion occurs as
      // described above, or until sqlite3_step() or sqlite3_reset() or
      // sqlite3_finalize() is called. The memory space used to hold strings
      // and BLOBs is freed automatically. Do not pass the pointers returned
      // from sqlite3_column_blob(), sqlite3_column_text(), etc. into
      // sqlite3_free().

      sqlite3_prepare(db, select_sql.c_str(), -1, &stmt, NULL);
      sqlite3_step( stmt );

      while (sqlite3_column_type( stmt, 0 ) != SQLITE_NULL)
      {
        MSChromatogram<> chrom;

        // int chrom_id = sqlite3_column_int(stmt, 0);
        const unsigned char * native_id = sqlite3_column_text(stmt, 1);
        chrom.setNativeID( std::string(reinterpret_cast<const char*>(native_id), sqlite3_column_bytes(stmt, 1)));
        String peptide_sequence;

        OpenMS::Precursor precursor;
        OpenMS::Product product;
        if (sqlite3_column_type(stmt, 2) != SQLITE_NULL) precursor.setCharge(sqlite3_column_int(stmt, 2));
        if (sqlite3_column_type(stmt, 3) != SQLITE_NULL) precursor.setDriftTime(sqlite3_column_double(stmt, 3));
        if (sqlite3_column_type(stmt, 4) != SQLITE_NULL) precursor.setMZ(sqlite3_column_double(stmt, 4));
        if (sqlite3_column_type(stmt, 5) != SQLITE_NULL) precursor.setIsolationWindowLowerOffset(sqlite3_column_double(stmt, 5));
        if (sqlite3_column_type(stmt, 6) != SQLITE_NULL) precursor.setIsolationWindowUpperOffset(sqlite3_column_double(stmt, 6));
        if (sqlite3_column_type(stmt, 7) != SQLITE_NULL) 
        {
          const unsigned char * pepseq = sqlite3_column_text(stmt, 7);
          peptide_sequence = std::string(reinterpret_cast<const char*>(pepseq), sqlite3_column_bytes(stmt, 7));
          precursor.setMetaValue("peptide_sequence", peptide_sequence);
        }
        // if (sqlite3_column_type(stmt, 8) != SQLITE_NULL) product.setCharge(sqlite3_column_int(stmt, 8));
        if (sqlite3_column_type(stmt, 9) != SQLITE_NULL) product.setMZ(sqlite3_column_double(stmt, 9));
        if (sqlite3_column_type(stmt, 10) != SQLITE_NULL) product.setIsolationWindowLowerOffset(sqlite3_column_double(stmt, 10));
        if (sqlite3_column_type(stmt, 11) != SQLITE_NULL) product.setIsolationWindowUpperOffset(sqlite3_column_double(stmt, 11));

        chrom.setPrecursor(precursor);
        chrom.setProduct(product);
        chromatograms.push_back(chrom);

        sqlite3_step( stmt );
      }

      // free memory
      sqlite3_finalize(stmt);
    }

};

class TOPPOpenSwathMzMLFileCacher
  : public TOPPBase,
    public ProgressLogger
{
 public:

  TOPPOpenSwathMzMLFileCacher()
    : TOPPBase("OpenSwathMzMLFileCacher","This tool caches the spectra and chromatogram data of an mzML to disk.", false)
  {
  }

 typedef MSExperiment<Peak1D> MapType;

 protected:

  void registerOptionsAndFlags_()
  {
    registerInputFile_("in","<file>","","Input mzML file");
    registerStringOption_("in_type", "<type>", "", "input file type -- default: determined from file extension or content\n", false);
    String formats("mzML,sqMass");
    setValidFormats_("in", ListUtils::create<String>(formats));
    setValidStrings_("in_type", ListUtils::create<String>(formats));

    formats = "mzML,sqMass";
    registerOutputFile_("out", "<file>", "", "Output file");
    setValidFormats_("out", ListUtils::create<String>(formats));
    registerStringOption_("out_type", "<type>", "", "Output file type -- default: determined from file extension or content\nNote: that not all conversion paths work or make sense.", false);
    setValidStrings_("out_type", ListUtils::create<String>(formats));

    //registerStringOption_("out_meta","<file>","","output file", false);
    //setValidFormats_("out_meta",ListUtils::create<String>("mzML"));

    registerFlag_("convert_back", "Convert back to mzML");
  }

  void convertToSQL(MSExperiment<> exp, String outputname)
  {
    SqMassWriter sql_mass(outputname);
    sql_mass.init();
    sql_mass.writeChromatograms(exp.getChromatograms());
    sql_mass.writeSpectra(exp.getSpectra());
  }

  void convertFromSQL(String inputname, MSExperiment<>& exp)
  {
    SqMassReader sql_mass_reader(inputname);
    sql_mass_reader.read(exp);
  }

  ExitCodes main_(int , const char**)
  {
    String out_meta = getStringOption_("out");
    String out_cached = out_meta + ".cached";
    bool convert_back =  getFlag_("convert_back");

    FileHandler fh;

    //input file type
    String in = getStringOption_("in");
    String in_cached = in + ".cached";
    FileTypes::Type in_type = FileTypes::nameToType(getStringOption_("in_type"));

    if (in_type == FileTypes::UNKNOWN)
    {
      in_type = fh.getType(in);
      writeDebug_(String("Input file type: ") + FileTypes::typeToName(in_type), 2);
    }

    if (in_type == FileTypes::UNKNOWN)
    {
      writeLog_("Error: Could not determine input file type!");
      return PARSE_ERROR;
    }

    //output file names and types
    String out = getStringOption_("out");
    FileTypes::Type out_type = FileTypes::nameToType(getStringOption_("out_type"));

    if (out_type == FileTypes::UNKNOWN)
    {
      out_type = fh.getTypeByFileName(out);
    }

    if (out_type == FileTypes::UNKNOWN)
    {
      writeLog_("Error: Could not determine output file type!");
      return PARSE_ERROR;
    }

    if (in_type == FileTypes::SQMASS && out_type == FileTypes::MZML)
    {
      MapType exp;
      convertFromSQL(in, exp);
      MzMLFile f;
      f.store(out, exp);
      return EXECUTION_OK;
    }
    else if (in_type == FileTypes::MZML && out_type == FileTypes::SQMASS)
    {
      std::cout << " out is sqmass, in is mzML" << std::endl;
      MzMLFile f;
      MapType exp;
      f.load(in, exp);
      convertToSQL(exp, out);
      return EXECUTION_OK;
    }


    if (!convert_back)
    {
      MapType exp;
      CachedmzML cacher;
      MzMLFile f;

      cacher.setLogType(log_type_);
      f.setLogType(log_type_);

      f.load(in,exp);
      cacher.writeMemdump(exp, out_cached);
      cacher.writeMetadata(exp, out_meta, true);
    }
    else
    {
      MzMLFile f;
      MapType meta_exp;
      CachedmzML cacher;
      MapType exp_reading;

      cacher.setLogType(log_type_);
      f.setLogType(log_type_);

      f.load(in,meta_exp);
      cacher.readMemdump(exp_reading, in_cached);

      std::cout << " read back, got " << exp_reading.size() << " spectra " << exp_reading.getChromatograms().size() << " chromats " << std::endl;

      {
      for (Size i=0; i<meta_exp.size(); ++i)
      {
        for (Size j = 0; j < meta_exp[i].getDataProcessing().size(); j++)
        {
          if (meta_exp[i].getDataProcessing()[j]->metaValueExists("cached_data"))
          {
            meta_exp[i].getDataProcessing()[j]->removeMetaValue("cached_data");
          }
        }
      }

      for (Size i=0; i < meta_exp.getNrChromatograms(); ++i)
      {
        for (Size j = 0; j < meta_exp.getChromatogram(i).getDataProcessing().size(); j++)
        {
          if (meta_exp.getChromatogram(i).getDataProcessing()[j]->metaValueExists("cached_data"))
          {
            meta_exp.getChromatogram(i).getDataProcessing()[j]->removeMetaValue("cached_data");
          }
        }
      }
      }

      if (meta_exp.size() != exp_reading.size())
      {
        std::cerr << " Both experiments need to have the same size!";
      }

      for (Size i=0; i<exp_reading.size(); ++i)
      {
        for (Size j = 0; j < exp_reading[i].size(); j++)
        {
          meta_exp[i].push_back(exp_reading[i][j]);
        }
      }
      std::vector<MSChromatogram<ChromatogramPeak> > chromatograms = exp_reading.getChromatograms();
      std::vector<MSChromatogram<ChromatogramPeak> > old_chromatograms = meta_exp.getChromatograms();
      for (Size i=0; i<chromatograms.size(); ++i)
      {
        for (Size j = 0; j < chromatograms[i].size(); j++)
        {
          old_chromatograms[i].push_back(chromatograms[i][j]);
        }
      }
      meta_exp.setChromatograms(old_chromatograms);


      f.store(out_meta,meta_exp);
    }

    return EXECUTION_OK;
  }
};

int main( int argc, const char** argv )
{

  TOPPOpenSwathMzMLFileCacher tool;
  return tool.main(argc,argv);
}

/// @endcond
