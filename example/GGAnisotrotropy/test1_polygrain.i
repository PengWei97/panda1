[Mesh]
  type = GeneratedMesh
  dim = 2 # Problem dimension
  nx = 10 # Number of elements in the x-direction
  ny = 10 # Number of elements in the y-direction
  nz = 0 # Number of elements in the z-direction
  xmin = 0    # minimum x-coordinate of the mesh
  xmax = 1000 # maximum x-coordinate of the mesh
  ymin = 0    # minimum y-coordinate of the mesh
  ymax = 1000 # maximum y-coordinate of the mesh
  zmin = 0
  zmax = 0
  elem_type = QUAD4 # Type of elements used in the mesh
  uniform_refine = 4 # Initial uniform refinement of the mesh

  parallel_type = replicated # Periodic BCs
[]

[GlobalParams]
  op_num = 10 # Number of grains
  var_name_base = gr # Base name of grains
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
    [./PolycrystalRandomIC]
      random_type = discrete
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
    T = 450 # Kelvin
    # molar_volume_value = 7.11e-6 # Units:m^3/mol
    Anisotropic_GB_file_name = anisotropy_mobility.txt   # anisotropy_energy.txt
    inclination_anisotropy = false # true
    # wGB = 10 # nm
    # length_scale = 1.0e-9
    # time_scale = 1.0e-9
    # delta_sigma = 0.1
    # delta_mob = 0.1
    wGB = 14 # Width of the diffuse GB
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
  [./num_grains]
    type = FeatureFloodCount
    variable = bnds
    threshold = 0.7
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
  end_time = 4000
  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.5 # Initial time step.  In this simulation it changes.
    optimal_iterations = 8 # Time step will adapt to maintain this number of nonlinear iterations
    growth_factor = 1.25
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
  file_base = test1_polygrain/polygrain
  execute_on = 'initial timestep_end'
  exodus = true
  csv = true
[]
