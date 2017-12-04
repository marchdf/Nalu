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
#include "TurbKineticEnergySSTDESSrcElemKernel.h"
#include "SpecificDissipationRateSSTSrcElemKernel.h"

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

namespace TurbKineticEnergySSTDESSrcElemKernel {

static constexpr double lhs[8][8] = {
  {
    0.18162464537165954, 0.060541548457219853, 0.020180516152406618,
    0.060541548457219853, 0.060541548457219853, 0.020180516152406618,
    0.0067268387174688722, 0.020180516152406618,
  },
  {
    0.045307794491460572, 0.13592338347438171, 0.045307794491460572,
    0.01510259816382019, 0.01510259816382019, 0.045307794491460572,
    0.01510259816382019, 0.0050341993879400634,
  },
  {
    0.012383713186335632, 0.037151139559006896, 0.11145341867702067,
    0.037151139559006896, 0.0041279043954452104, 0.012383713186335632,
    0.037151139559006896, 0.012383713186335632,
  },
  {
    0.050131830062481619, 0.016710610020827209, 0.050131830062481619,
    0.15039549018744486, 0.016710610020827209, 0.0055702033402757357,
    0.016710610020827209, 0.050131830062481619,
  },
  {
    0.045307794491460572, 0.01510259816382019, 0.0050341993879400634,
    0.01510259816382019, 0.13592338347438171, 0.045307794491460572,
    0.01510259816382019, 0.045307794491460572,
  },
  {
    0.01090404787216135, 0.03271214361648405, 0.01090404787216135,
    0.0036346826240537832, 0.03271214361648405, 0.098136430849452144,
    0.03271214361648405, 0.01090404787216135,
  },
  {
    0.0029833026770613251, 0.008949908031183975, 0.026849724093551925,
    0.008949908031183975, 0.008949908031183975, 0.026849724093551925,
    0.080549172280655779, 0.026849724093551925,
  },
  {
    0.012383713186335632, 0.0041279043954452104, 0.012383713186335632,
    0.037151139559006896, 0.037151139559006896, 0.012383713186335632,
    0.037151139559006896, 0.11145341867702067,
  },
};

static constexpr double rhs[8] = {
  -0.90508598874801927, -0.63112775201881044, -0.56515879124551993,
  -0.83337912275421244, -0.65875459678226767, -0.43406521038693352,
  -0.37461216408191811, -0.57407140728772554,
};

} // namespace TurbKineticEnergySSTDESSrcElemKernel

namespace SpecificDissipationRateSSTSrcElemKernel {

static constexpr double lhs[8][8] = {
  {
    0.16830144483151507, 0.056100481610505029, 0.018700160536835007,
    0.056100481610505029, 0.056100481610505029, 0.018700160536835007,
    0.0062333868456116697, 0.018700160536835007,
  },
  {
    0.039427962222187668, 0.118283886666563, 0.039427962222187668,
    0.013142654074062557, 0.013142654074062557, 0.039427962222187668,
    0.013142654074062557, 0.0043808846913541855,
  },
  {
    0.011393016582828068, 0.034179049748484201, 0.10253714924545261,
    0.034179049748484201, 0.0037976721942760226, 0.011393016582828068,
    0.034179049748484201, 0.011393016582828068,
  },
  {
    0.048595301433219079, 0.016198433811073026, 0.048595301433219079,
    0.14578590429965724, 0.016198433811073026, 0.0053994779370243424,
    0.016198433811073026, 0.048595301433219079,
  },
  {
    0.038838931011900821, 0.012946310337300274, 0.0043154367791000915,
    0.012946310337300274, 0.11651679303570248, 0.038838931011900821,
    0.012946310337300274, 0.038838931011900821,
  },
  {
    0.0075958254763728835, 0.022787476429118651, 0.0075958254763728835,
    0.0025319418254576278, 0.022787476429118651, 0.068362429287355947,
    0.022787476429118651, 0.0075958254763728835,
  },
  {
    0.0022113769883583762, 0.0066341309650751286, 0.019902392895225385,
    0.0066341309650751286, 0.0066341309650751286, 0.019902392895225385,
    0.059707178685676154, 0.019902392895225385,
  },
  {
    0.01041806031306038, 0.0034726867710201267, 0.01041806031306038,
    0.031254180939181142, 0.031254180939181142, 0.01041806031306038,
    0.031254180939181142, 0.093762542817543426,
  },
};

static constexpr double rhs[8] = {
  0.74767435661938209, 0.54339860453045752, 0.48953955267320354,
  0.69657914067554871, 0.52854536355285653, 0.31844406143891979,
  0.30053506738269653, 0.48071925627889456,
};

} // namespace SpecificDissipationRateSSTSrcElemKernel
} // namespace hex8_golds
} // anonymous namespace

TEST_F(SSTKernelHex8Mesh, turbkineticenergysstsrcelem)
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

TEST_F(SSTKernelHex8Mesh, turbkineticenergysstdessrcelem)
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
    new sierra::nalu::TurbKineticEnergySSTDESSrcElemKernel<
      sierra::nalu::AlgTraitsHex8>(
      bulk_, solnOpts_, helperObjs.assembleElemSolverAlg->dataNeededByKernels_,
      false));

  // Add to kernels to be tested
  helperObjs.assembleElemSolverAlg->activeKernels_.push_back(kernel.get());

  helperObjs.assembleElemSolverAlg->execute();

  EXPECT_EQ(helperObjs.linsys->lhs_.dimension(0), 8u);
  EXPECT_EQ(helperObjs.linsys->lhs_.dimension(1), 8u);
  EXPECT_EQ(helperObjs.linsys->rhs_.dimension(0), 8u);

  namespace gold_values = hex8_golds::TurbKineticEnergySSTDESSrcElemKernel;
  unit_test_kernel_utils::expect_all_near(
    helperObjs.linsys->rhs_, gold_values::rhs);
  unit_test_kernel_utils::expect_all_near<8>(
    helperObjs.linsys->lhs_, gold_values::lhs);
}

TEST_F(SSTKernelHex8Mesh, specificdissipationratesstsrcelem)
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
    new sierra::nalu::SpecificDissipationRateSSTSrcElemKernel<
      sierra::nalu::AlgTraitsHex8>(
      bulk_, solnOpts_, helperObjs.assembleElemSolverAlg->dataNeededByKernels_,
      false));

  // Add to kernels to be tested
  helperObjs.assembleElemSolverAlg->activeKernels_.push_back(kernel.get());

  helperObjs.assembleElemSolverAlg->execute();

  EXPECT_EQ(helperObjs.linsys->lhs_.dimension(0), 8u);
  EXPECT_EQ(helperObjs.linsys->lhs_.dimension(1), 8u);
  EXPECT_EQ(helperObjs.linsys->rhs_.dimension(0), 8u);

  namespace gold_values = hex8_golds::SpecificDissipationRateSSTSrcElemKernel;
  unit_test_kernel_utils::expect_all_near(
    helperObjs.linsys->rhs_, gold_values::rhs);
  unit_test_kernel_utils::expect_all_near<8>(
    helperObjs.linsys->lhs_, gold_values::lhs);
}
