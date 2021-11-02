# This example tests the implementation of PolycrstalStoredEnergy kernels that assigns excess stored energy to grains with dislocation density
[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 96
  ny = 32
  nz = 0
  xmin = 0
  xmax = 96
  ymin = 0
  ymax = 32
[]

[GlobalParams]
  block = 0
  op_num = 12
  deformed_grain_num = 25
  var_name_base = gr
  grain_num = 30
  grain_tracker = grain_tracker
  time_scale = 1e-2
  length_scale = 1e-8
[]

[Variables]
  [./PolycrystalVariables]
  [../]
[]

[AuxVariables]
  [./bnds]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[UserObjects]
  [./voronoi]
    type = PolycrystalVoronoi
    rand_seed = 81
    coloring_algorithm = jp
  [../]
  [./grain_tracker]
    type = GrainTracker
    threshold = 0.2
    connecting_threshold = 0.08
    compute_var_to_feature_map = true
    flood_entity_type = elemental
    execute_on = ' initial timestep_begin'
    outputs = none
  [../]
[]

[ICs]
  [./PolycrystalICs]
    [./PolycrystalColoringIC]
      polycrystal_ic_uo = voronoi
    [../]
  [../]
[]

[BCs]
  [./Periodic]
    [./all]
      auto_direction = 'x y'
    [../]
  [../]
[]

[Kernels]
  [./PolycrystalKernel]
  [../]
  [./PolycrystalStoredEnergy]
    # ACSEDGPoly -- adds the stored energy contribution to grain growth models.
    # op_num = N, 创建N个ACSEDGPoly kernels
    # var_name_base = gr
    # deformed_grain_num

  [../]
[]

[AuxKernels]
  [./BndsCalc]
    type = BndsCalcAux
    variable = bnds
    execute_on = timestep_end
  [../]
[]

[Materials]
  [./deformed]
    type = DeformedGrainMaterial
    int_width = 4.0
    outputs = exodus
    # outputs: kappa_op, gamma_asymm,..., Disloc_Den_i(disloc density), 
  [../]
[]

[Preconditioning]
  [./SMP]
    type = SMP
    full = true
  [../]
[]

[Executioner]
  type = Transient
  nl_max_its = 15
  scheme = bdf2
  solve_type = PJFNK
  petsc_options_iname = -pc_type
  petsc_options_value = asm
  l_max_its = 15
  l_tol = 1.0e-3
  nl_rel_tol = 1.0e-8
  start_time = 0.0
  end_time = 200
  nl_abs_tol = 1e-8
  # dt = 0.20
  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.20
    growth_factor = 1.2
    cutback_factor = 0.8
    optimal_iterations = 8
  [../]
  [./Adaptivity]
    initial_adaptivity = 2
    refine_fraction = 0.8
    coarsen_fraction = 0.05
    max_h_level = 3
  [../]
[]

[Outputs]
  file_base = ./deformation/out_deformated
  exodus = true
  interval = 1
  show = bnds
  perf_graph = true
[]
