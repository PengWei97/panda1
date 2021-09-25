[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 10
  ny = 3
  xmax = 1000
  ymax = 1000
  elem_type = QUAD4
  uniform_refine = 2
[]

[GlobalParams]
  op_num = 2
  var_name_base = gr
[]

[Variables]
  [./PolycrystalVariables]
     order = FIRST
     family = LAGRANGE
  [../]
[]

[ICs]
  [./PolycrystalICs]
    [./BicrystalBoundingBoxIC]
      x1 = 0
      y1 = 0
      x2 = 500
      y2 = 1000
    [../]
  [../]
[]

[AuxVariables]
  [./bnds]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[Kernels]
  [./PolycrystalKernel]
  [../]
[]

[AuxKernels]
  [./bnds_aux]
    type = BndsCalcAux
    variable = bnds
    execute_on = timestep_end
  [../]
[]

[BCs]
  # [./top_displacement]
  #   type = DirichletBC
  #   variable = disp_y
  #   boundary = top
  #   value = -10.0
  # [../]
  # [./x_anchor]
  #   type = DirichletBC
  #   variable = disp_x
  #   boundary = 'left right'
  #   value = 0.0
  # [../]
  # [./y_anchor]
  #   type = DirichletBC
  #   variable = disp_y
  #   boundary = bottom
  #   value = 0.0
  # [../]
[]

[Materials]
  [./Copper]
    type = GBEvolution
    block = 0
    T = 500 # K
    wGB = 75 # nm
    GBmob0 = 2.5e-6 #m^4/(Js) from Schoenfelder 1997
    Q = 0.23 #Migration energy in eV
    GBenergy = 0.708 #GB energy in J/m^2
    time_scale = 1.0e-6
    length_scale = 1.0e-9
  [../]
  # [./ElasticityTensor]
  #   type = ComputePolycrystalElasticityTensor
  #   grain_tracker = grain_tracker
  # [../]
  # [./strain]
  #   type = ComputeSmallStrain
  #   block = 0
  #   displacements = 'disp_x disp_y'
  # [../]
  # [./stress]
  #   type = ComputeLinearElasticStress
  #   block = 0
  # [../]
[]

[UserObjects]
  # [./euler_angle_file]
  #   type = EulerAngleFileReader
  #   file_name = test.tex
  # [../]
  # [./grain_tracker]
  #   type = GrainTrackerElasticity
  #   connecting_threshold = 0.05
  #   compute_var_to_feature_map = true
  #   flood_entity_type = elemental
  #   execute_on = 'initial timestep_begin'

  #   euler_angle_provider = euler_angle_file
  #   fill_method = symmetric9
  #   C_ijkl = '1.27e5 0.708e5 0.708e5 1.27e5 0.708e5 1.27e5 0.7355e5 0.7355e5 0.7355e5'

  #   outputs = none
  # [../]
[]

[Postprocessors]
  [./dt]
    type = TimestepSize
  [../]
  [./gr0_area]
    type = ElementIntegralVariablePostprocessor
    variable = gr0
  [../]
[]

[Preconditioning]
  [./SMP]
   type = SMP
   coupled_groups = 'gr0,gr1'
  [../]
[]

[Executioner]
  type = Transient

  solve_type = 'PJFNK'
  petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart -pc_hypre_boomeramg_strong_threshold'
  petsc_options_value = 'hypre boomeramg 31 0.7'

  l_max_its = 30
  l_tol = 1e-4
  nl_max_its = 30
  nl_rel_tol = 1e-9

  start_time = 0.0
  num_steps = 3
  dt = 0.2

  [./Adaptivity]
   initial_adaptivity = 2
    refine_fraction = 0.7
    coarsen_fraction = 0.1
    max_h_level = 2
  [../]
[]

[Outputs]
  execute_on = 'timestep_end'
  exodus = true
[]

