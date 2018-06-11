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

#include <OpenMS/ANALYSIS/OPENSWATH/DIAHelper.h>

#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/CoarseIsotopePatternGenerator.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderAlgorithmPickedHelperStructs.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderAlgorithm.h>

#include <OpenMS/MATH/MISC/CubicSpline2d.h>
#include <OpenMS/MATH/MISC/SplineBisection.h>
#include <OpenMS/FILTERING/SMOOTHING/GaussFilterAlgorithm.h>

#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/MSExperiment.h>

#include <iostream>
#include <utility>
#include <boost/bind.hpp>

namespace OpenMS
{
  namespace DIAHelpers
  {


    void fitSplineToPeak(OpenSwath::SpectrumPtr spectrum, const double left, const double right,
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
      std::vector<double>::iterator central_mz_it;
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
