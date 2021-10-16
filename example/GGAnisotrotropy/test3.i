[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 40
  ny = 40
  nz = 0
  xmin = 0
  xmax = 1000
  ymin = 0
  ymax = 1000
  zmin = 0
  zmax = 0
  elem_type = QUAD4
[]

[GlobalParams]
  op_num = 2
  var_name_base = gr
  wGB = 10
  length_scale = 1.0e-9
  time_scale = 1.0e-9
[]

[Variables]
  [./PolycrystalVariables]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[ICs]
  [./PolycrystalICs]
    [./BicrystalCircleGrainIC]
      radius = 333.33
      x = 500
      y = 500
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
    Anisotropic_GB_file_name = anisotropy.txt

    inclination_anisotropy = true
    delta_sigma = 0.1
    delta_mob = 0.0
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

  start_time = 0.0
  end_time = 300.0
  # dt = 10.0
  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 10 # Initial time step.  In this simulation it changes.
    optimal_iterations = 6 # Time step will adapt to maintain this number of nonlinear iterations
  [../]

  [./Adaptivity]
    # Block that turns on mesh adaptivity. Note that mesh will never coarsen beyond initial mesh (before uniform refinement)
    initial_adaptivity = 3 # Number of times mesh is adapted to initial condition
    refine_fraction = 0.7 # Fraction of high error that will be refined
    coarsen_fraction = 0.1 # Fraction of low error that will coarsened
    max_h_level = 5 # Max number of refinements used, starting from initial mesh (before uniform refinement)
  [../]
[]

[Outputs]
  file_base = test3_out/test3
  execute_on = 'timestep_end'
  exodus = true
  csv = true
[]
