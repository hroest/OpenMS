// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2021.
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

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/DATASTRUCTURES/Param.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/CHEMISTRY/ResidueModification.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <iostream>
#include <omp.h>

using namespace OpenMS;
using namespace std;

Int main()
{

  cout << String("Tutorial is working") << endl;

  auto p = Param();
  cout << p << endl;

   static ModificationsDB* mdb = ModificationsDB::getInstance();

   int nr_iterations (1e4), test (0);
	 omp_set_num_threads(8);
   // std::cout << "Setting up nested loop with " << omp_get_max_threads() << " threads " << std::endl;
#pragma omp parallel for reduction (+: test)
  for (int k = 1; k < nr_iterations + 1; k++)
  {
    int mod_id = k;
    String modname = "mod" + String(mod_id);
    std::unique_ptr<ResidueModification> new_mod(new ResidueModification());
    new_mod->setFullId(modname);
    new_mod->setMonoMass( 0.11 * mod_id);
    new_mod->setAverageMass(1.0);
    new_mod->setDiffMonoMass( 0.05 * mod_id);
      mdb->addModification(std::move(new_mod));
			int tmp = (int)mdb->getModification(modname)->getAverageMass();
    	test += tmp;
		  // std::cout << "tmp " << tmp << " with thread " << omp_get_thread_num() << " it " << k<< std::endl;
  }

  std::cout << "got a total of " << test << std::endl;

  AASequence seq = AASequence::fromString("PEPTIDE");

  // throw Exception::NotImplemented(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);

  return 0;
} //end of main

