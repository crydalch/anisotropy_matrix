/* ************************************************************************
 * Copyright 2013 Alexander Mishurov
 * 
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 * http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 * ***********************************************************************
*/

#ifndef HDKPLUGIN_SOPANISOTROPICMATRIX_H_
#define HDKPLUGIN_SOPANISOTROPICMATRIX_H_

#include <omp.h>
#include <SYS/SYS_Math.h>
#include <UT/UT_DSOVersion.h>
#include <UT/UT_Interrupt.h>
#include <UT/UT_Matrix4.h>
#include <OP/OP_Operator.h>
#include <OP/OP_OperatorTable.h>
#include <SOP/SOP_Node.h>
#include <PRM/PRM_Include.h>
#include <GU/GU_Detail.h>
#include <GEO/GEO_PointTree.h>
#include "GEO/GEO_PrimPoly.h"
#include "GU/GU_PrimSphere.h"
#include "UT/UT_Matrix3.h"

#include <GU/GU_PrimPacked.h>

// for threading
#include <UT/UT_ParallelUtil.h>
#include <GA/GA_PageIterator.h>
#include <GA/GA_PageHandle.h>
#include <UT/UT_LockUtil.h>
#include <UT/UT_Lock.h>

namespace HDK_AMPlugins {

class SOP_AnisotropyMatrix : public SOP_Node
{
    public:
      SOP_AnisotropyMatrix(OP_Network *, const char *, OP_Operator *);
      virtual ~SOP_AnisotropyMatrix();
      static OP_Node *myConstructor(OP_Network *, const char *, OP_Operator *);
    public:
      static PRM_Template myTemplateList[];

    protected:
      virtual OP_ERROR cookMySop(OP_Context &);
};

class covarianceMatrixTask {
    public:
        covarianceMatrixTask(   GU_Detail *myParticleGdp, GU_Detail *newSphereGdp,/* GA_Attribute *attr_particle_p,*/
                                fpreal &smoothing_kernel_radius,
                                fpreal &search_radius,
                                fpreal &scale_addition,
                                unsigned particles_threshold,
                                UT_AutoInterrupt &boss,
                                unsigned write_attr_only):

        myParticleGdp(myParticleGdp),
        newSphereGdp(newSphereGdp),
        /*attr_particle_p(attr_particle_p),*/
        smoothing_kernel_radius(smoothing_kernel_radius),
        search_radius(search_radius),
        scale_addition(scale_addition),
        particles_threshold(particles_threshold),
        boss(boss),
        write_attr_only(write_attr_only)
        
    {}

    void operator()(const GA_SplittableRange &sr) const;

    private:
        GU_Detail *myParticleGdp;
        GU_Detail *newSphereGdp;
        /*GA_Attribute *attr_particle_p;*/
        fpreal &smoothing_kernel_radius;
        fpreal &search_radius;
        fpreal &scale_addition;
        unsigned particles_threshold;
        UT_AutoInterrupt &boss;
        unsigned write_attr_only;
};

} // HDK_AMPlugins namespace

#endif  // HDKPLUGIN_SOPANISOTROPICMATRIX_H_