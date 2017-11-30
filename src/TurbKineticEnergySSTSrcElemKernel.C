/*------------------------------------------------------------------------*/
/*  Copyright 2014 National Renewable Energy Laboratory.                  */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#include "TurbKineticEnergySSTSrcElemKernel.h"
#include "FieldTypeDef.h"
#include "SolutionOptions.h"

#include "BuildTemplates.h"
#include "ScratchViews.h"

// stk_mesh/base/fem
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>

namespace sierra {
  namespace nalu {

    template<typename AlgTraits>
    TurbKineticEnergySSTSrcElemKernel<AlgTraits>::TurbKineticEnergySSTSrcElemKernel(
										    const stk::mesh::BulkData& bulkData,
										    const SolutionOptions& solnOpts,
										    ElemDataRequests& dataPreReqs)
      : Kernel(),
	betaStar_(solnOpts.get_turb_model_constant(TM_betaStar)),
	tkeProdLimitRatio_(solnOpts.get_turb_model_constant(TM_tkeProdLimitRatio)),
	ipNodeMap_(sierra::nalu::MasterElementRepo::get_volume_master_element(AlgTraits::topo_)->ipNodeMap())
    {
      const stk::mesh::MetaData& metaData = bulkData.mesh_meta_data();
      ScalarFieldType *tke = metaData.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "turbulent_ke");
      tkeNp1_ = &tke->field_of_state(stk::mesh::StateNP1);
      ScalarFieldType *sdr = metaData.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "specific_dissipation_rate");
      sdrNp1_ = &sdr->field_of_state(stk::mesh::StateNP1);
      ScalarFieldType *density = metaData.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "density");
      densityNp1_ = &density->field_of_state(stk::mesh::StateNP1);
      VectorFieldType *velocity = metaData.get_field<VectorFieldType>(stk::topology::NODE_RANK, "velocity");
      velocityNp1_ = &(velocity->field_of_state(stk::mesh::StateNP1));
      tvisc_ = metaData.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "turbulent_viscosity");
      coordinates_ = metaData.get_field<VectorFieldType>(stk::topology::NODE_RANK, solnOpts.get_coordinates_name());

      MasterElement* meSCV = sierra::nalu::MasterElementRepo::get_volume_master_element(AlgTraits::topo_);

      // add master elements
      dataPreReqs.add_cvfem_volume_me(meSCV);

      // fields and data
      dataPreReqs.add_coordinates_field(*coordinates_, AlgTraits::nDim_, CURRENT_COORDINATES);
      dataPreReqs.add_gathered_nodal_field(*tkeNp1_, 1);
      dataPreReqs.add_gathered_nodal_field(*sdrNp1_, 1);
      dataPreReqs.add_gathered_nodal_field(*densityNp1_, 1);
      dataPreReqs.add_gathered_nodal_field(*velocityNp1_, AlgTraits::nDim_);
      dataPreReqs.add_gathered_nodal_field(*tvisc_, 1);
      dataPreReqs.add_master_element_call(SCV_VOLUME, CURRENT_COORDINATES);
      dataPreReqs.add_master_element_call(SCV_GRAD_OP, CURRENT_COORDINATES);
    }

    template<typename AlgTraits>
    TurbKineticEnergySSTSrcElemKernel<AlgTraits>::~TurbKineticEnergySSTSrcElemKernel()
    {}

    template<typename AlgTraits>
    void
    TurbKineticEnergySSTSrcElemKernel<AlgTraits>::execute(
							  SharedMemView<DoubleType**>& lhs,
							  SharedMemView<DoubleType *>& rhs,
							  ScratchViews<DoubleType>& scratchViews)
    {

      SharedMemView<DoubleType*>& v_tkeNp1 = scratchViews.get_scratch_view_1D(*tkeNp1_);
      SharedMemView<DoubleType*>& v_sdrNp1 = scratchViews.get_scratch_view_1D(*sdrNp1_);
      SharedMemView<DoubleType*>& v_densityNp1 = scratchViews.get_scratch_view_1D(*densityNp1_);
      SharedMemView<DoubleType**>& v_velocityNp1 = scratchViews.get_scratch_view_2D(*velocityNp1_);
      SharedMemView<DoubleType*>& v_tvisc = scratchViews.get_scratch_view_1D(*tvisc_);
      SharedMemView<DoubleType***>& v_dndx = scratchViews.get_me_views(CURRENT_COORDINATES).dndx_scv;
      SharedMemView<DoubleType*>& v_scv_volume = scratchViews.get_me_views(CURRENT_COORDINATES).scv_volume;

      for (int ip=0; ip < AlgTraits::numScvIp_; ++ip) {

	// nearest node to ip
	const int nearestNode = ipNodeMap_[ip];

	// save off scvol
	const DoubleType scV = v_scv_volume(ip);

	for (int ic=0; ic < AlgTraits::nodesPerElement_; ++ic) {

	  const DoubleType rhoIC = v_densityNp1(ic);
	  const DoubleType sdrIC = v_sdrNp1(ic);
	  const DoubleType tkeIC = v_tkeNp1(ic);
	  const DoubleType tviscIC = v_tvisc(ic);

	  DoubleType Pk = 0.0;
	  for ( int i = 0; i < AlgTraits::nDim_; ++i ) {
	    const DoubleType dni = v_dndx(ip,ic,i);
	    const DoubleType ui = v_velocityNp1(ic,i);
	    for ( int j = 0; j < AlgTraits::nDim_; ++j ) {
	      const DoubleType dnj = v_dndx(ip,ic,j);
	      Pk += dnj * ui * (dnj*ui + dni * v_velocityNp1(ic,j));
	   }
	  }
	  Pk *= tviscIC;

	  // tke factor
	  const DoubleType tkeFac = betaStar_*rhoIC*sdrIC;

	  // dissipation and production (limited)
	  DoubleType Dk = tkeFac * tkeIC;
	  Pk = stk::math::min(Pk, tkeProdLimitRatio_*Dk);
	  
	  // assemble RHS and LHS
	  rhs(nearestNode) += (Pk - Dk)*scV;
	  lhs(nearestNode,nearestNode) += tkeFac*scV;
	}
      }
    }

    INSTANTIATE_KERNEL(TurbKineticEnergySSTSrcElemKernel);

  } // nalu
} // sierra
