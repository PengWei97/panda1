[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 12
  ny = 12
  nz = 0
  xmin = 0
  xmax = 1000
  ymin = 0
  ymax = 1000
  zmin = 0
  zmax = 0
  uniform_refine = 3
  elem_type = QUAD4
[]

[GlobalParams]
  op_num = 15
  var_name_base = gr
  wGB = 14
  length_scale = 1.0e-9
  time_scale = 1.0e-9
[]

[UserObjects]
  [./voronoi]
    type = PolycrystalVoronoi
    grain_num = 15
    rand_seed = 42
    coloring_algorithm = bt # We must use bt to force the UserObject to assign one grain to each op
  [../]
[]

[Variables]
  [./PolycrystalVariables]
  [../]
[]

[ICs]
  [./PolycrystalICs]
    [./PolycrystalColoringIC]
      polycrystal_ic_uo = voronoi
    [../]
  [../]
[]

[AuxVariables]
  [./bnds]
    order = FIRST
    family = LAGRANGE
  [../]
  [./unique_grains]
    order = FIRST
    family = LAGRANGE
  [../]
  [./var_indices]
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
  [./Periodic]
    [./top_bottom]
      auto_direction = 'x y'
    [../]
  [../]
[]

[Materials]
  [./CuGrGranisotropic]
    type = GBAnisotropy
    T = 600 # K

    # molar_volume_value = 7.11e-6 #Units:m^3/mol
    Anisotropic_GB_file_name = poly_anisotropy_mobility.txt   # anisotropy_energy.txt
    inclination_anisotropy = false # true
    outputs = exodus
  [../]
[]

[Postprocessors]
  [./dt]
    # Outputs the current time step
    type = TimestepSize
  [../]

  [./gr1_area]
    type = ElementIntegralVariablePostprocessor
    variable = gr1
  [../]
  [./gr2_area]
    type = ElementIntegralVariablePostprocessor
    variable = gr2
  [../]
[]

[Executioner]
  type = Transient
  scheme = bdf2

  #Preconditioned JFNK (default)
  solve_type = 'PJFNK'

  petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart'
  petsc_options_value = 'hypre boomeramg 31'

  l_max_its = 30
  l_tol = 1e-4
  nl_max_its = 40
  nl_rel_tol = 1e-9

  start_time = 0
  end_time = 1000
  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 5
    growth_factor = 1.2
    cutback_factor = 0.8
    optimal_iterations = 8
  [../]
  [./Adaptivity]
    initial_adaptivity = 4
    refine_fraction = 0.8
    coarsen_fraction = 0.05
    max_h_level = 3
  [../]
[]

[Outputs]
  file_base = ./poly_test1/poly_test1_GBEnergy
  execute_on = 'timestep_end'
  interval = 3
  exodus = true
  csv = true
[]
