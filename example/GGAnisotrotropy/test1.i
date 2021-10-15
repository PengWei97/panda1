[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 60
  ny = 30
  nz = 0
  xmin = 0
  xmax = 1000
  ymin = 0
  ymax = 600
  zmin = 0
  zmax = 0
  elem_type = QUAD4
[]

[GlobalParams]
  op_num = 3
  var_name_base = gr
  wGB = 10
  length_scale = 1.0e-9
  time_scale = 1.0e-9
[]

[Variables]
  [./PolycrystalVariables]
    # action
    # op_num = 3
    # var_name_base = gr
    # gr0 gr1 gr3
  [../]
[]

[ICs]
  [./PolycrystalICs]
    [./Tricrystal2CircleGrainsIC]
      # action
      # op_num = 3 
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
    # molar_volume_value = 7.11e-6 # Units:m^3/mol
    Anisotropic_GB_file_name = anisotropy_mobility.txt   # anisotropy_energy.txt
    inclination_anisotropy = false # true
    # wGB = 50 
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

  num_steps = 20
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
  file_base = outputs/test1
  execute_on = 'timestep_end'
  exodus = true
  csv = true
[]
