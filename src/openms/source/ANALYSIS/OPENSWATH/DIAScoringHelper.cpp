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
// $Maintainer: Witold Wolski, Hannes Roest $
// $Authors: Witold Wolski, Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/DIAScoringHelper.h>

#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/CoarseIsotopePatternGenerator.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderAlgorithmPickedHelperStructs.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderAlgorithm.h>

#include <OpenMS/MATH/MISC/CubicSpline2d.h>
#include <OpenMS/MATH/MISC/SplineBisection.h>
#include <OpenMS/FILTERING/SMOOTHING/GaussFilterAlgorithm.h>

#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/CONCEPT/Macros.h>

#include <iostream>
#include <utility>
#include <boost/bind.hpp>

namespace OpenMS
{
  namespace DIAHelpers
  {

    void adjustExtractionWindow(double& right, double& left, const double& mz_extraction_window, const bool& dia_extraction_ppm_)
    {
      OPENMS_PRECONDITION(mz_extraction_window > 0, "MZ extraction window needst to be larger than zero.");

      if (dia_extraction_ppm_)
      {
        left -= left * mz_extraction_window / 2e6;
        right += right * mz_extraction_window / 2e6;
      }
      else
      {
        left -= mz_extraction_window / 2.0;
        right += mz_extraction_window / 2.0;
      }
    }

    void fitSplineToPeak(OpenSwath::SpectrumPtr /* spectrum */, const double /* left */, const double /* right */,
                         const std::vector<double> & newmz, const std::vector<double> & newint,
                         double& max_peak_mz )
    {
#if 0
      std::vector<double> fnewmz;
      std::vector<double> fnewint;
      GaussFilterAlgorithm f;
      fnewmz.resize(newmz.size());
      fnewint.resize(newint.size());
      // void initialize(double gaussian_width, double spacing, double ppm_tolerance, bool use_ppm_tolerance);
      f.initialize(10, 0.001, 10, true); // TODO algorithm params!
      // f.initialize(10, 0.01, 1, true); // TODO algorithm params!
      f.initialize(10, 0.001, 20, true); // TODO algorithm params!
      f.filter< std::vector<double>::const_iterator, std::vector<double>::iterator >(newmz.begin(), newmz.end(), newint.begin(), fnewmz.begin(), fnewint.begin());
#else
      std::vector<double> fnewmz = newmz;
      std::vector<double> fnewint = newint;
#endif

      std::map<double, double> peak_raw_data;
      size_t maxk = -1;

      for (Size k = 0; k < fnewmz.size(); k++)
      {
        peak_raw_data[ fnewmz[k] ] = fnewint[k];
        if (fnewint[maxk] < fnewint[k])
        {
          maxk = k;
        }

      }

      // std::cout << " peak reaw data " << peak_raw_data.size() << std::endl;
      max_peak_mz = -1;
      if (peak_raw_data.size() < 3) {return;}
      CubicSpline2d peak_spline (peak_raw_data);

      // calculate maximum by evaluating the spline's 1st derivative
      // (bisection method)
      max_peak_mz = fnewmz[maxk];
      double max_peak_int = fnewint[maxk];
      double threshold = 1e-6;
      double left_neighbor_mz = fnewmz[ std::max( (int)maxk-1, 0)];
      double right_neighbor_mz = fnewmz[std::min( (Size)fnewmz.size()-1, maxk+1)];
      OpenMS::Math::spline_bisection(peak_spline, left_neighbor_mz, right_neighbor_mz, max_peak_mz, max_peak_int, threshold);
    }

    void integrateWindows(const OpenSwath::SpectrumPtr spectrum,
                          const std::vector<double> & windowsCenter, double width,
                          std::vector<double> & integratedWindowsIntensity,
                          std::vector<double> & integratedWindowsMZ,
                          bool remZero)
    {
      std::vector<double>::const_iterator beg = windowsCenter.begin();
      std::vector<double>::const_iterator end = windowsCenter.end();
      double mz, intensity;
      for (; beg != end; ++beg)
      {
        double left = *beg - width / 2.0;
        double right = *beg + width / 2.0;
        if (integrateWindow(spectrum, left, right, mz, intensity, false))
        {
          integratedWindowsIntensity.push_back(intensity);
          integratedWindowsMZ.push_back(mz);
        }
        else if (!remZero)
        {
          integratedWindowsIntensity.push_back(0.);
          integratedWindowsMZ.push_back(*beg);
        }
        else
        {
        }
      }
    }

    void integrateDriftSpectrum(OpenSwath::SpectrumPtr spectrum, 
                                double mz_start,
                                double mz_end,
                                double & im,
                                double & intensity,
                                std::vector<std::pair<double, double> >& res, 
                                double drift_start,
                                double drift_end)
    {
      OPENMS_PRECONDITION(spectrum->getDriftTimeArray() != nullptr, "Cannot filter by drift time if no drift time is available.");
      OPENMS_PRECONDITION(spectrum->getMZArray()->data.size() == spectrum->getIntensityArray()->data.size(), "MZ and Intensity array need to have the same length.");
      OPENMS_PRECONDITION(spectrum->getMZArray()->data.size() == spectrum->getDriftTimeArray()->data.size(), "MZ and Drift Time array need to have the same length.");
      OPENMS_PRECONDITION(std::adjacent_find(spectrum->getMZArray()->data.begin(),
              spectrum->getMZArray()->data.end(), std::greater<double>()) == spectrum->getMZArray()->data.end(),
              "Precondition violated: m/z vector needs to be sorted!" )

      im = 0;
      intensity = 0;

      // rounding multiplier for the ion mobility value
      // TODO: how to improve this -- will work up to 42949.67296
      double IM_IDX_MULT = 10e5;

      std::map< int, double> im_chrom;
      {
        // get the weighted average for noncentroided data.
        // TODO this is not optimal if there are two peaks in this window (e.g. if the window is too large)
        typedef std::vector<double>::const_iterator itType;

        itType mz_arr_end = spectrum->getMZArray()->data.end();
        itType int_it = spectrum->getIntensityArray()->data.begin();
        itType im_it = spectrum->getDriftTimeArray()->data.begin();

        // this assumes that the spectra are sorted!
        itType mz_it = std::lower_bound(spectrum->getMZArray()->data.begin(),
          spectrum->getMZArray()->data.end(), mz_start);
        itType mz_it_end = std::lower_bound(mz_it, mz_arr_end, mz_end);

        // also advance intensity and ion mobility iterator now
        std::iterator_traits< itType >::difference_type iterator_pos = std::distance((itType)spectrum->getMZArray()->data.begin(), mz_it);
        std::advance(int_it, iterator_pos);
        std::advance(im_it, iterator_pos);

        // Iterate from mz start to end, only storing ion mobility values that are in the range
        for (; mz_it != mz_it_end; ++mz_it, ++int_it, ++im_it)
        {
          if ( *im_it >= drift_start && *im_it <= drift_end)
          {
            // std::cout << "IM " << *im_it << " mz " << *mz_it << " int " << *int_it << std::endl;
            im_chrom[ int((*im_it)*IM_IDX_MULT) ] += *int_it;
            intensity += (*int_it);
            im += (*int_it) * (*im_it);
          }
        }

        if (intensity > 0.)
        {
          // std::cout << " before " << im << std::endl;
          im /= intensity;
          // std::cout << " after " << im << std::endl;
        }
        else
        {
          im = -1;
          intensity = 0;
        }

      }

      for (auto k : im_chrom) 
      {
        res.push_back(std::make_pair( k.first / IM_IDX_MULT, k.second ) );
      }

    }

    void integrateDriftSpectrum_x(OpenSwath::SpectrumPtr spectrum, 
                                              double mz_start,
                                              double mz_end,
                                              double & im,
                                              double & intensity,
                                              double drift_start,
                                              double drift_end)
    {
      OPENMS_PRECONDITION(spectrum->getDriftTimeArray() != nullptr, "Cannot filter by drift time if no drift time is available.");
      OPENMS_PRECONDITION(spectrum->getMZArray()->data.size() == spectrum->getIntensityArray()->data.size(), "MZ and Intensity array need to have the same length.");
      OPENMS_PRECONDITION(spectrum->getMZArray()->data.size() == spectrum->getDriftTimeArray()->data.size(), "MZ and Drift Time array need to have the same length.");
      OPENMS_PRECONDITION(std::adjacent_find(spectrum->getMZArray()->data.begin(),
              spectrum->getMZArray()->data.end(), std::greater<double>()) == spectrum->getMZArray()->data.end(),
              "Precondition violated: m/z vector needs to be sorted!" )

      im = 0;
      intensity = 0;

      // get the weighted average for noncentroided data.
      // TODO this is not optimal if there are two peaks in this window (e.g. if the window is too large)
      typedef std::vector<double>::const_iterator itType;

      itType mz_arr_end = spectrum->getMZArray()->data.end();
      itType int_it = spectrum->getIntensityArray()->data.begin();
      itType im_it = spectrum->getDriftTimeArray()->data.begin();

      // this assumes that the spectra are sorted!
      itType mz_it = std::lower_bound(spectrum->getMZArray()->data.begin(),
          spectrum->getMZArray()->data.end(), mz_start);
      itType mz_it_end = std::lower_bound(mz_it, mz_arr_end, mz_end);

      // also advance intensity and ion mobility iterator now
      std::iterator_traits< itType >::difference_type iterator_pos = std::distance((itType)spectrum->getMZArray()->data.begin(), mz_it);
      std::advance(int_it, iterator_pos);
      std::advance(im_it, iterator_pos);

      // Iterate from mz start to end, only storing ion mobility values that are in the range
      for (; mz_it != mz_it_end; ++mz_it, ++int_it, ++im_it)
      {
        if ( *im_it >= drift_start && *im_it <= drift_end)
        {
          intensity += (*int_it);
          im += (*int_it) * (*im_it);
        }
      }

      if (intensity > 0.)
      {
        im /= intensity;
      }
      else
      {
        im = -1;
        intensity = 0;
      }

    }

    bool integrateWindow(const OpenSwath::SpectrumPtr spectrum,
                         double mz_start,
                         double mz_end,
                         double & mz,
                         double & intensity,
                         bool centroided)
    {
      OPENMS_PRECONDITION(spectrum->getMZArray()->data.size() == spectrum->getIntensityArray()->data.size(), "MZ and Intensity array need to have the same length.");
      OPENMS_PRECONDITION(std::adjacent_find(spectrum->getMZArray()->data.begin(),
              spectrum->getMZArray()->data.end(), std::greater<double>()) == spectrum->getMZArray()->data.end(),
              "Precondition violated: m/z vector needs to be sorted!" )

      mz = 0;
      intensity = 0;

      if (!centroided)
      {
        // get the weighted average for noncentroided data.
        // TODO this is not optimal if there are two peaks in this window (e.g. if the window is too large)
        typedef std::vector<double>::const_iterator itType;

        itType mz_arr_end = spectrum->getMZArray()->data.end();
        itType int_it = spectrum->getIntensityArray()->data.begin();

        // this assumes that the spectra are sorted!
        itType mz_it = std::lower_bound(spectrum->getMZArray()->data.begin(),
          spectrum->getMZArray()->data.end(), mz_start);
        itType mz_it_end = std::lower_bound(mz_it, mz_arr_end, mz_end);

        // also advance intensity iterator now
        std::iterator_traits< itType >::difference_type iterator_pos = std::distance((itType)spectrum->getMZArray()->data.begin(), mz_it);
        std::advance(int_it, iterator_pos);

        for (; mz_it != mz_it_end; ++mz_it, ++int_it)
        {
          intensity += (*int_it);
          mz += (*int_it) * (*mz_it);
        }

        if (intensity > 0.)
        {
          mz /= intensity;
          return true;
        }
        else
        {
          mz = -1;
          intensity = 0;
          return false;
        }

      }
      else
      {
        // not implemented
        throw "Not implemented";
      }
    }


    // for SWATH -- get the theoretical b and y series masses for a sequence
    void getBYSeries(const AASequence& a, //
                     std::vector<double>& bseries, //
                     std::vector<double>& yseries, //
                     TheoreticalSpectrumGenerator const * generator,
                     UInt charge)
    {
      OPENMS_PRECONDITION(charge > 0, "Charge is a positive integer");
      //too slow!
      //TheoreticalSpectrumGenerator generator;
      //Param p;
      //p.setValue("add_metainfo", "true",
      //           "Adds the type of peaks as metainfo to the peaks, like y8+, [M-H2O+2H]++");
      //generator.setParameters(p);
      PeakSpectrum spec;
      generator->getSpectrum(spec, a, charge, charge);

      const PeakSpectrum::StringDataArray& ion_name = spec.getStringDataArrays()[0];

      for (Size i = 0; i != spec.size(); ++i)
      {
        if (ion_name[i][0] == 'y')
        {
          yseries.push_back(spec[i].getMZ());
        }
        else if (ion_name[i][0] == 'b')
        {
          bseries.push_back(spec[i].getMZ());
        }
      }
    } // end getBYSeries

    // for SWATH -- get the theoretical b and y series masses for a sequence
    void getTheorMasses(const AASequence& a, std::vector<double>& masses,
                        TheoreticalSpectrumGenerator const * generator,
                        UInt charge)
    {
      OPENMS_PRECONDITION(charge > 0, "Charge is a positive integer");
      //too slow!
      //TheoreticalSpectrumGenerator generator;
      //Param p;
      //p.setValue("add_metainfo", "false",
      //           "Adds the type of peaks as metainfo to the peaks, like y8+, [M-H2O+2H]++");
      //p.setValue("add_precursor_peaks", "true", "Adds peaks of the precursor to the spectrum, which happen to occur sometimes");
      //generator.setParameters(p);
      PeakSpectrum spec;
      generator->getSpectrum(spec, a, charge, charge);
      for (PeakSpectrum::iterator it = spec.begin();
           it != spec.end(); ++it)
      {
        masses.push_back(it->getMZ());
      }
    } // end getBYSeries

    void getAveragineIsotopeDistribution(const double product_mz,
                                         std::vector<std::pair<double, double> >& isotopesSpec, const double charge,
                                         const int nr_isotopes, const double mannmass)
    {
      typedef OpenMS::FeatureFinderAlgorithmPickedHelperStructs::TheoreticalIsotopePattern TheoreticalIsotopePattern;
      // create the theoretical distribution
      CoarseIsotopePatternGenerator solver(nr_isotopes);
      TheoreticalIsotopePattern isotopes;
      //std::cout << product_mz * charge << std::endl;
      auto d = solver.estimateFromPeptideWeight(product_mz * charge);

      double mass = product_mz;
      for (IsotopeDistribution::Iterator it = d.begin(); it != d.end(); ++it)
      {
        isotopesSpec.push_back(std::make_pair(mass, it->getIntensity()));
        mass += mannmass;
      }
    } //end of dia_isotope_corr_sub

    //simulate spectrum from AASequence
    void simulateSpectrumFromAASequence(const AASequence& aa,
                                        std::vector<double>& firstIsotopeMasses, //[out]
                                        std::vector<std::pair<double, double> >& isotopeMasses, //[out]
                                        TheoreticalSpectrumGenerator const * generator, double charge)
    {
      getTheorMasses(aa, firstIsotopeMasses, generator, charge);
      for (std::size_t i = 0; i < firstIsotopeMasses.size(); ++i)
      {
        getAveragineIsotopeDistribution(firstIsotopeMasses[i], isotopeMasses,
                                        charge);
      }
    }

    //given an experimental spectrum add isotope pattern.
    void addIsotopes2Spec(const std::vector<std::pair<double, double> >& spec,
                          std::vector<std::pair<double, double> >& isotopeMasses, //[out]
                          double charge)
    {

      for (std::size_t i = 0; i < spec.size(); ++i)
      {
        std::vector<std::pair<double, double> > isotopes;
        getAveragineIsotopeDistribution(spec[i].first, isotopes, charge);
        for (Size j = 0; j < isotopes.size(); ++j)
        {
          isotopes[j].second *= spec[i].second; //multiple isotope intensity by spec intensity
          isotopeMasses.push_back(isotopes[j]);
        }
      }
    }

    //Add masses before first isotope
    void addPreisotopeWeights(const std::vector<double>& firstIsotopeMasses,
                              std::vector<std::pair<double, double> >& isotopeSpec, // output
                              UInt nrpeaks, double preIsotopePeaksWeight, // weight of pre isotope peaks
                              double mannmass, double charge)
    {
      for (std::size_t i = 0; i < firstIsotopeMasses.size(); ++i)
      {
        double mul = 1.;
        for (UInt j = 0; j < nrpeaks; ++j, ++mul)
        {
          isotopeSpec.push_back(
            std::make_pair(firstIsotopeMasses[i] - (mul * mannmass) / charge,
                           preIsotopePeaksWeight));
        }
      }
      sortByFirst(isotopeSpec);
    }

    struct MassSorter :
      std::binary_function<double, double, bool>
    {
      bool operator()(const std::pair<double, double>& left,
                      const std::pair<double, double>& right)
      {
        return left.first < right.first;
      }
    };

    void sortByFirst(std::vector<std::pair<double, double> >& tmp)
    {
      std::sort(tmp.begin(), tmp.end(), MassSorter());
    }

    void extractFirst(const std::vector<std::pair<double, double> >& peaks,
                      std::vector<double>& mass)
    {
      for (std::size_t i = 0; i < peaks.size(); ++i)
      {
        mass.push_back(peaks[i].first);
      }
    }

    void extractSecond(const std::vector<std::pair<double, double> >& peaks,
                       std::vector<double>& mass)
    {
      for (std::size_t i = 0; i < peaks.size(); ++i)
      {
        mass.push_back(peaks[i].second);
      }
    }

    //modify masses by charge
    void modifyMassesByCharge(
      const std::vector<std::pair<double, double> >& isotopeSpec,
      std::vector<std::pair<double, double> >& resisotopeSpec, double charge)
    {
      resisotopeSpec.clear();
      std::pair<double, double> tmp_;
      for (std::size_t i = 0; i < isotopeSpec.size(); ++i)
      {
        tmp_ = isotopeSpec[i];
        tmp_.first /= charge;
        resisotopeSpec.push_back(tmp_);
      }
    }

//simulate spectrum from AASequence
//void simulateSpectrumFromTransitions(AASequence & aa,
//std::vector<double> & firstIsotopeMasses,
//std::vector<std::pair<double, double> > & isotopeMasses, uint32_t charge)

  }
}
