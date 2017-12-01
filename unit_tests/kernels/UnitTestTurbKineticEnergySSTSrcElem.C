/*------------------------------------------------------------------------*/
/*  Copyright 2014 National Renewable Energy Laboratory.                  */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#include "kernels/UnitTestKernelUtils.h"
#include "UnitTestUtils.h"
#include "UnitTestHelperObjects.h"

#include "TurbKineticEnergySSTSrcElemKernel.h"

TEST_F(TurbKineticEnergySSTKernelHex8Mesh, turbkineticenergysstsrcelem)
{

  fill_mesh_and_init_fields();

  // Setup solution options
  solnOpts_.meshMotion_ = false;
  solnOpts_.meshDeformation_ = false;
  solnOpts_.externalMeshDeformation_ = false;
  solnOpts_.initialize_turbulence_constants();

  unit_test_utils::HelperObjectsNewME helperObjs(bulk_, stk::topology::HEX_8, 1, partVec_[0]);


  // Initialize the kernel
  std::unique_ptr<sierra::nalu::Kernel> kernel(
					       new sierra::nalu::TurbKineticEnergySSTSrcElemKernel<sierra::nalu::AlgTraitsHex8>(
																bulk_, solnOpts_, helperObjs.assembleElemSolverAlg->dataNeededByKernels_, false)
  );

  // Add to kernels to be tested
  helperObjs.assembleElemSolverAlg->activeKernels_.push_back(kernel.get());

  helperObjs.assembleElemSolverAlg->execute();

  EXPECT_EQ(helperObjs.linsys->lhs_.dimension(0), 8u);
  EXPECT_EQ(helperObjs.linsys->lhs_.dimension(1), 8u);
  EXPECT_EQ(helperObjs.linsys->rhs_.dimension(0), 8u);

}
