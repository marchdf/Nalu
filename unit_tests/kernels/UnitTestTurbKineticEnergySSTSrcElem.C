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

namespace {
namespace hex8_golds {
namespace TurbKineticEnergySSTSrcElemKernel {

static constexpr double lhs[8][8] = {
  {
    0.0069752671206081166, 0.0023250890402027055, 0.00077502968006756855,
    0.0023250890402027055, 0.0023250890402027055, 0.00077502968006756855,
    0.00025834322668918952, 0.00077502968006756855,
  },
  {
    0.0017833984706739155, 0.005350195412021746, 0.0017833984706739155,
    0.00059446615689130516, 0.00059446615689130516, 0.0017833984706739155,
    0.00059446615689130516, 0.00019815538563043504,
  },
  {
    0.00047062403152685068, 0.0014118720945805519, 0.0042356162837416562,
    0.0014118720945805519, 0.00015687467717561689, 0.00047062403152685068,
    0.0014118720945805519, 0.00047062403152685068,
  },
  {
    0.00185533509705029, 0.00061844503235009669, 0.00185533509705029,
    0.0055660052911508696, 0.00061844503235009669, 0.00020614834411669889,
    0.00061844503235009669, 0.00185533509705029,
  },
  {
    0.0018553350970502902, 0.0006184450323500968, 0.00020614834411669891,
    0.0006184450323500968, 0.0055660052911508705, 0.0018553350970502902,
    0.0006184450323500968, 0.0018553350970502902,
  },
  {
    0.00047062403152685068, 0.0014118720945805519, 0.00047062403152685068,
    0.00015687467717561689, 0.0014118720945805519, 0.0042356162837416562,
    0.0014118720945805519, 0.00047062403152685068,
  },
  {
    0.00013065390263746531, 0.00039196170791239594, 0.0011758851237371878,
    0.00039196170791239594, 0.00039196170791239594, 0.0011758851237371878,
    0.0035276553712115634, 0.0011758851237371878,
  },
  {
    0.00052603049446954737, 0.00017534349815651577, 0.00052603049446954737,
    0.0015780914834086419, 0.0015780914834086419, 0.00052603049446954737,
    0.0015780914834086419, 0.0047342744502259261,
  },
};

static constexpr double rhs[8] = {
  -0.0097567470924261324, 0.026675107608750099,   0.038686221727213671,
  0.020801003295001737,   -0.0020389462344089851, 0.032579280131220176,
  0.043395843349824381,   0.026965189169389121,
};

} // namespace TurbKineticEnergySSTSrcElemKernel
} // namespace hex8_golds
} // anonymous namespace

TEST_F(TurbKineticEnergySSTKernelHex8Mesh, turbkineticenergysstsrcelem)
{

  fill_mesh_and_init_fields();

  // Setup solution options
  solnOpts_.meshMotion_ = false;
  solnOpts_.meshDeformation_ = false;
  solnOpts_.externalMeshDeformation_ = false;
  solnOpts_.initialize_turbulence_constants();

  unit_test_utils::HelperObjectsNewME helperObjs(
    bulk_, stk::topology::HEX_8, 1, partVec_[0]);

  // Initialize the kernel
  std::unique_ptr<sierra::nalu::Kernel> kernel(
    new sierra::nalu::TurbKineticEnergySSTSrcElemKernel<
      sierra::nalu::AlgTraitsHex8>(
      bulk_, solnOpts_, helperObjs.assembleElemSolverAlg->dataNeededByKernels_,
      false));

  // Add to kernels to be tested
  helperObjs.assembleElemSolverAlg->activeKernels_.push_back(kernel.get());

  helperObjs.assembleElemSolverAlg->execute();

  EXPECT_EQ(helperObjs.linsys->lhs_.dimension(0), 8u);
  EXPECT_EQ(helperObjs.linsys->lhs_.dimension(1), 8u);
  EXPECT_EQ(helperObjs.linsys->rhs_.dimension(0), 8u);

  namespace gold_values = hex8_golds::TurbKineticEnergySSTSrcElemKernel;
  unit_test_kernel_utils::expect_all_near(
    helperObjs.linsys->rhs_, gold_values::rhs);
  unit_test_kernel_utils::expect_all_near<8>(
    helperObjs.linsys->lhs_, gold_values::lhs);
}
