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
 * ************************************************************************/

#include "SOP_anisotropy_matrix.h"

using namespace HDK_AMPlugins;

void
newSopOperator(OP_OperatorTable *table)
{
  table->addOperator(new OP_Operator("anisotropy_matrix",
                                     "AnisotropyMatrix",
                                     SOP_AnisotropyMatrix::myConstructor,
                                     SOP_AnisotropyMatrix::myTemplateList,
                                     1,
                                     1,
                                     0));
}

static PRM_Name krnl_name("kernel", "Kernel Radius");
static PRM_Name srch_name("search", "Search Radius");
static PRM_Name scl_name("scale", "Scale Addition");
static PRM_Name thrs_name("threshold", "Particles Threshold");
static PRM_Name usethreads_name("usethreads", "Threaded");
static PRM_Name attronly_name("attronly", "Point Attr Only");

static PRM_Default krnl_default(1);
static PRM_Default srch_default(2);
static PRM_Default scl_default(0.1);
static PRM_Default thrs_default(6);
static PRM_Default usethreads_default(1);
static PRM_Default attronly_default(1);

static PRM_Range unsigned_range(PRM_RANGE_RESTRICTED, 0);
static PRM_Range unsigned_10_range(PRM_RANGE_UI, 0, PRM_RANGE_RESTRICTED, 10);
static PRM_Range unsigned_50_range(PRM_RANGE_UI, 0, PRM_RANGE_RESTRICTED, 50);

PRM_Template
SOP_AnisotropyMatrix::myTemplateList[] = {
    PRM_Template(PRM_TOGGLE,  1, &usethreads_name, &attronly_default),
    PRM_Template(PRM_TOGGLE,  1, &attronly_name, &usethreads_default),
    PRM_Template(PRM_FLT,  1, &krnl_name, &krnl_default, 0, &unsigned_10_range),
    PRM_Template(PRM_FLT,  1, &srch_name, &srch_default, 0, &unsigned_10_range),
    PRM_Template(PRM_FLT,  1, &scl_name, &scl_default, 0, &unsigned_range),
    PRM_Template(PRM_INT,  1, &thrs_name, &thrs_default, 0, &unsigned_50_range),
    PRM_Template(),
};

OP_Node *
SOP_AnisotropyMatrix::myConstructor(OP_Network *net,
                                    const char *name,
                                    OP_Operator *op)
{
    return new SOP_AnisotropyMatrix(net, name, op);
}

SOP_AnisotropyMatrix::SOP_AnisotropyMatrix(OP_Network *net,
                                           const char *name,
                                           OP_Operator *op)
: SOP_Node(net, name, op) {}

SOP_AnisotropyMatrix::~SOP_AnisotropyMatrix() {}


void
covarianceMatrixTask::
operator()(const GA_SplittableRange &sr) const
{
    // Attribute handles for P and aniso_mtx
    GA_ROPageHandleV3   hndl_geo_p(myParticleGdp->getP()); // Why does this work, but not passing the GA_Attribute?
    GA_RWPageHandleM3D  hndl_geo_amtx(myParticleGdp->addFloatTuple(GA_ATTRIB_POINT, "aniso_mtx", 9));

    // Create a KD Tree for nearest neighbour lookups
    GEO_PointTreeGAOffset tree;
    tree.build(myParticleGdp);//->getPos3AsArray());

    // Generate a UT_Lock, to ensure our new GU_Detail is only written to by one thread at a time
    static UT_Lock theMutex;

    // Iterate over the pages allocated to a given thread
    for (GA_PageIterator pit = sr.beginPages(); !pit.atEnd(); ++pit)
    {
        if (!boss.wasInterrupted())
        {
            // Iterate over the blocks of each page
            GA_Offset block_offset_start, block_offset_end;
            for(GA_Iterator it(pit.begin()); it.blockAdvance(block_offset_start, block_offset_end);)
            {
                if (!boss.wasInterrupted())
                {
                    // Iterate over the offsets of each block
                    hndl_geo_p.setPage(block_offset_start);
                    hndl_geo_amtx.setPage(block_offset_start);
                    for (GA_Offset ptoff = block_offset_start; ptoff < block_offset_end; ++ptoff)
                    {
                        //UT_Vector3 particle_pos = hndl_geo_p.get(ptoff);
                        UT_Vector3  particle_pos = myParticleGdp->getPos3(ptoff);

                        // Find the points within search_radius, save them to an array
                        GEO_PointTreeGAOffset::IdxArrayType close_particles_indices;
                        tree.findAllCloseIdx(   particle_pos, 
                                                search_radius,
                                                close_particles_indices);

                        unsigned close_particles_count = close_particles_indices.entries();

                        UT_Matrix3 anisotropy_matrix;

                        if (close_particles_count > 0) 
                        {   

                            // Calculation of weighted mean
                            UT_Vector3 weighted_mean(0, 0, 0);
                            fpreal weight = 0;
                            fpreal weighting_function = 0;
                            UT_Vector3 weighted_position(0, 0, 0);

                            for (unsigned i = 0; i < close_particles_count; i++) 
                            {
                                //UT_Vector3 close_particle_pos = hndl_geo_p.get(close_particles_indices(i));
                                UT_Vector3 close_particle_pos = myParticleGdp->getPos3(close_particles_indices(i));

                                UT_Vector3 distance = particle_pos - close_particle_pos;
                                weight = 1 - std::pow((distance.length() / search_radius), 3);

                                weighting_function += weight;
                                weighted_position += close_particle_pos * weight;
                            }

                            if (weighting_function != 0)
                                weighted_mean = weighted_position/weighting_function;
                            else
                                weighted_mean = weighted_position;

                            // Calculation of covariance matrix and SVD -- example code provided by ndickson, thank you!
                            UT_Matrix3D covariance_matrix(0); 

                            for (unsigned i = 0; i < close_particles_count; i++) 
                            { 
                                //UT_Vector3 close_particle_pos = hndl_geo_p.get(close_particles_indices(i));
                                UT_Vector3 close_particle_pos = myParticleGdp->getPos3(close_particles_indices(i));

                                UT_Vector3 weighted_distance = close_particle_pos - weighted_mean; 
                                UT_Vector3 distance = particle_pos - close_particle_pos; 
                                weight = 1 - std::pow((distance.length() / search_radius), 3); 
                                UT_Vector3 weighted = weight * weighted_distance; 

                                // Only 6 unique components, since symmetric 
                                covariance_matrix(0,0) += weighted(0)*weighted_distance(0); 
                                covariance_matrix(0,1) += weighted(0)*weighted_distance(1); 
                                covariance_matrix(0,2) += weighted(0)*weighted_distance(2); 
                                covariance_matrix(1,1) += weighted(1)*weighted_distance(1); 
                                covariance_matrix(1,2) += weighted(1)*weighted_distance(2); 
                                covariance_matrix(2,2) += weighted(2)*weighted_distance(2); 
                            } 

                            // Copy symmetric components 
                            covariance_matrix(1,0) = covariance_matrix(0,1); 
                            covariance_matrix(2,0) = covariance_matrix(0,2); 
                            covariance_matrix(2,1) = covariance_matrix(1,2); 

                            if (weighting_function != 0) 
                                covariance_matrix /= weighting_function; 

                            UT_Matrix3D rotation_matrix; 
                            UT_Matrix3D diagonal_matrix; // diagonal will hold eigenvalues 

                            // This is probably overkill; there are likely faster algorithms, but this should work, up to some tolerance. 
                            covariance_matrix.diagonalizeSymmetric(rotation_matrix,diagonal_matrix); 

                            UT_Vector3D eigen_values_vector(diagonal_matrix(0,0), diagonal_matrix(1,1), diagonal_matrix(2,2));
                            // End SVD/CovarianceMatrix

                            // Particles threshold
                            if (close_particles_count > particles_threshold) {
                                for (unsigned i = 0; i < 3; i++)
                                    eigen_values_vector(i) = 2 * eigen_values_vector(i) + scale_addition;
                            } 
                            else {
                                eigen_values_vector(0) =
                                eigen_values_vector(1) =
                                eigen_values_vector(2) = 1;
                            }

                            // Convert to HDK matrices
                            UT_Matrix3 rotation_matrix_hdk(rotation_matrix(0, 0),
                                                           rotation_matrix(0, 1),
                                                           rotation_matrix(0, 2),
                                                           rotation_matrix(1, 0),
                                                           rotation_matrix(1, 1),
                                                           rotation_matrix(1, 2),
                                                           rotation_matrix(2, 0),
                                                           rotation_matrix(2, 1),
                                                           rotation_matrix(2, 2));

                            UT_Matrix3 eigen_values_matrix_hdk(eigen_values_vector(0), 0, 0,
                                                               0, eigen_values_vector(1), 0,
                                                               0, 0, eigen_values_vector(2));

                            rotation_matrix_hdk.invert();
                            UT_Matrix3 rotation_matrix_transpose_hdk = rotation_matrix_hdk;
                            rotation_matrix_transpose_hdk.transpose();

                            // Compose anisotropy matrix
                            anisotropy_matrix = rotation_matrix_hdk *
                                                eigen_values_matrix_hdk *
                                                rotation_matrix_transpose_hdk;

                            anisotropy_matrix *= 1 / smoothing_kernel_radius;
                        }

                        // If we're not writing the attribute only, generate the spheres
                        if (!write_attr_only)
                        {
                            GU_Detail   *geometry_gdp = new GU_Detail;
                            GA_Offset tptoff = geometry_gdp->appendPoint();
                            GU_PrimSphereParms parms(geometry_gdp, tptoff);
                            parms.freq = 2;
                            parms.type = GEO_PATCH_TRIANGLE;
                            parms.xform.scale(0.2,0.2,0.2);
                            (GEO_PrimSphere*) GU_PrimSphere::build(parms, GEO_PRIMPOLY);

                            // Apply transform for each point to a sphere
                            GA_Offset       optoff;
                            GA_FOR_ALL_PTOFF(geometry_gdp, optoff)
                            {
                                UT_Vector3 geometry_pos = geometry_gdp->getPos3(optoff);
                                if (close_particles_count > 0)
                                    geometry_pos.colVecMult(anisotropy_matrix);
                                geometry_gdp->setPos3(optoff,geometry_pos + particle_pos);
                               
                            }

                            // Add geometry copy to final geometry
                            UT_Lock::Scope lock(theMutex);

                            newSphereGdp->copy( *geometry_gdp,
                                                GEO_COPY_ADD,
                                                true,
                                                true,
                                                GA_DATA_ID_BUMP);

                            UT_Lock::Scope unlock(theMutex);        
                        }
                        else
                            hndl_geo_amtx.set(ptoff, anisotropy_matrix);

                    }
                }//if (!boss.wasInterrupted()
            }
        }//if (!boss.wasInterrupted())
    }
};


OP_ERROR
SOP_AnisotropyMatrix::cookMySop(OP_Context &context)
{
    fpreal now = context.getTime();

    if (lockInputs(context) >= UT_ERROR_ABORT)
        return error();

    setupLocalVars();

    if (error() < UT_ERROR_ABORT) 
    {

        //UT_AutoInterrupt *progress = new UT_AutoInterrupt("Calculating matrices");
        UT_AutoInterrupt progress("Calculating matrices");

        GU_Detail   *particles_gdp_orig = (GU_Detail *)inputGeo(0, context);
        GU_Detail   *new_gdp = new GU_Detail;
        //GU_Detail *geometry_gdp = (GU_Detail *)inputGeo(1, context);

        //GU_Detail *temporary_gdp = (GU_Detail *)inputGeo(0, context);
        // GU_Detail pack_gdp;

        gdp->clearAndDestroy();
        //temporary_gdp->clearAndDestroy();

        // MATRIX PARMS
        fpreal smoothing_kernel_radius_pval = evalFloat("kernel", 0, now);
        fpreal search_radius_pval = evalFloat("search", 0, now);
        fpreal scale_addition_pval = evalFloat("scale", 0, now);
        unsigned particles_threshold_pval = evalInt("threshold", 0, now);
        unsigned use_threads = evalInt("usethreads", 0, now);
        unsigned attr_only = evalInt("attronly", 0, now);

        //GA_Attribute *attr_geo_p = particles_gdp_orig->getP();
        
        // Create the aniso_mtx attribute
        particles_gdp_orig->addFloatTuple(GA_ATTRIB_POINT, "aniso_mtx", 9);

        ///Try building inside the threaded kernel?
        //GEO_PointTreeGAOffset tree;
        //tree.build(particles_gdp_orig);//->getPos3AsArray());

        // for (unsigned i=0; i < particles_gdp_orig->getNumPoints(); i++)
        // {
        //     GA_Offset poffset = pack_gdp.appendPointOffset();
        //     GU_PrimPacked::build(pack_gdp, "PackedSphere", poffset);   
        // }

        const GA_SplittableRange sr(particles_gdp_orig->getPointRange());
        
        covarianceMatrixTask threadedCMtask(particles_gdp_orig, new_gdp, /*attr_geo_p,*/
                                            smoothing_kernel_radius_pval,
                                            search_radius_pval,
                                            scale_addition_pval,
                                            particles_threshold_pval, progress,
                                            attr_only);

        //covarianceMatrixTask threadedCMtask(particles_gdp_orig);
        if (use_threads)
            UTparallelForLightItems(sr,threadedCMtask);
        else
            UTserialFor(sr,threadedCMtask);

        if (attr_only)
        {
            gdp->copy( *particles_gdp_orig,
                        GEO_COPY_ADD,
                        true,
                        true,
                        GA_DATA_ID_BUMP);
        }
        else
        {
            gdp->copy( *new_gdp,
                        GEO_COPY_ADD,
                        true,
                        true,
                        GA_DATA_ID_BUMP);
        }

    }
    
    unlockInputs();
    resetLocalVarRefs();

    return error();
}
