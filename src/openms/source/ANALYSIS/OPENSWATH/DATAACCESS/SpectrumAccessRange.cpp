// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
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

#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SpectrumAccessRange.h>

#include <iostream>

namespace OpenMS
{

  SpectrumAccessRange::SpectrumAccessRange(
      OpenSwath::SpectrumAccessPtr sptr,
      double mz_low, double mz_high) :
        SpectrumAccessTransforming(sptr), 
        mz_low_(mz_low), 
        mz_high_(mz_high)
    {}
        
    SpectrumAccessRange::~SpectrumAccessRange() {}

    boost::shared_ptr<OpenSwath::ISpectrumAccess> SpectrumAccessRange::lightClone() const
    {
      // Create a light clone of *this by initializing a new
      // SpectrumAccessRange with a light clone of the underlying
      // SpectrumAccess object and the parameters.
      return boost::shared_ptr<SpectrumAccessRange>(
          new SpectrumAccessRange(sptr_->lightClone(), mz_low_, mz_high_));
    }

    OpenSwath::BinaryDataArrayPtr deepCopy(OpenSwath::BinaryDataArrayPtr d_)
    {
      OpenSwath::BinaryDataArrayPtr d (new OpenSwath::BinaryDataArray);
      d->data = d_->data;
      d->description = d_->description;
      return d;
    }

    OpenSwath::SpectrumPtr deepCopy(OpenSwath::SpectrumPtr s_)
    {
      std::vector<OpenSwath::BinaryDataArrayPtr> tmp;
      for (auto const &p: s_->getDataArrays())
      {
        tmp.push_back(deepCopy(p));
      }

      OpenSwath::SpectrumPtr s (new OpenSwath::Spectrum);
      s->getDataArrays() = tmp;
      return s;
    }

    OpenSwath::SpectrumPtr SpectrumAccessRange::getSpectrumById(int id)
    {
      // NOTE: we may have gotten a reference - we should not change the
      // underlying data structure but rather modify a copy of it.
      OpenSwath::SpectrumPtr s = deepCopy(sptr_->getSpectrumById(id));

      // create filter array which stores whether to remove the datapoint
      // (default = true = remove it)
      std::vector<int> tmp(s->getMZArray()->data.size(), true);
      auto s_it = tmp.begin();
      for (const auto& it : s->getMZArray()->data)
      {
        if ( it >= mz_low_ && it <= mz_high_) {*s_it = false;}
        s_it++;
      }

      // now filter all data arrays
      for (auto& arr : s->getDataArrays())
      {
        s_it = tmp.begin();
        arr->data.erase(
            std::remove_if(arr->data.begin(), arr->data.end(),
              [&s_it](const double &)
              {
                // return current iterator value and increment *after*
                return (*(s_it++));
              }),
            arr->data.end());
      }

      return s;
    }
}
