基于ykvishal，https://github.com/idaholab/moose/discussions/17422
# [Mesh]
#   type = GeneratedMesh
#   dim = 3
#   nx = 45
#   ny = 45
#   nz = 45
#   xmin = 0
#   xmax = 45
#   ymin = 0
#   ymax = 45
#   zmin = 0
#   zmax = 45
#   elem_type = HEX8
#   # parallel_type = distributed
# []

[Mesh]
  [gmg]
    type = DistributedRectilinearMeshGenerator
    dim = 2 #3
    nx = 2500 #180
    ny = 2500 #180
    #nz = 180
    xmin = 0
    xmax = 2500 #180
    ymin = 0
    ymax = 2500 #180
    #zmin = 0
    #zmax = 180
    elem_type = QUAD4  #HEX8
  []
[]


[GlobalParams]
  op_num = 28
  var_name_base = gr
[]

[Variables]
  [./PolycrystalVariables]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[UserObjects]
  [./voronoi]
    type = PolycrystalVoronoi
    grain_num = 100 # 6000 # Number of grains
    rand_seed = 200 # 301
    coloring_algorithm = jp
    int_width = 6 # Width of diffuse interfaces
    # 0 for 2 days, 6 for 24 hrs
    # if Projecting initial condition takes 50% more time
  [../]
  [./term]
    type = Terminator
    expression = 'grain_tracker < 10' # 218
  [../]
  [./grain_tracker]
    type = GrainTracker
    threshold = 0.2
    connecting_threshold = 0.08
    compute_halo_maps = true # Only necessary for displaying HALOS
    compute_var_to_feature_map = true
    polycrystal_ic_uo = voronoi
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
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./ghost_elements]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./halos]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./var_indices]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[Kernels]
  [./PolycrystalKernel]
  [../]
[]

[AuxKernels]
  [./BndsCalc]
    type = BndsCalcAux
    variable = bnds
    execute_on = 'initial timestep_end'
  [../]
  [./unique_grains]
    type = FeatureFloodCountAux
    variable = unique_grains
    field_display = UNIQUE_REGION
    execute_on = 'initial timestep_end'
    flood_counter = grain_tracker
  [../]
  [./ghost_elements]
    type = FeatureFloodCountAux
    variable = ghost_elements
    field_display = GHOSTED_ENTITIES
    execute_on = 'initial timestep_end'
    flood_counter = grain_tracker
  [../]
  [./halos]
    type = FeatureFloodCountAux
    variable = halos
    field_display = HALOS
    execute_on = 'initial timestep_end'
    flood_counter = grain_tracker
  [../]
  [./var_indices]
    type = FeatureFloodCountAux
    variable = var_indices
    field_display = VARIABLE_COLORING
    execute_on = 'initial timestep_end'
    flood_counter = grain_tracker
  [../]
[]

[BCs]
 [./Periodic]
   [./All]
     auto_direction = 'x y'
   [../]
 [../]
[]

[Materials]
  [./Copper]
    type = GBEvolution
    T = 500
    wGB = 3 # um
    GBmob0 = 2.5e-6 #m^4/(Js) from Schoenfelder 1997
    Q = 0.23 #Migration energy in eV
    GBenergy = 0.708 #GB energy in J/m^2
    molar_volume = 7.11e-6 #Molar volume in m^3/mol
    length_scale = 1.0e-6
    time_scale = 1.0
  [../]
[]

[Postprocessors]
  [./dt]
    type = TimestepSize
  [../]
  [./n_elements]
    type = NumElems
    execute_on = timestep_end
  [../]
  [./n_nodes]
    type = NumNodes
    execute_on = timestep_end
  [../]
  [./DOFs]
    type = NumDOFs
  [../]
  # [./grain_tracker]
  #   type = GrainTracker
  #   threshold = 0.1
  #   compute_halo_maps = true
  # [../]
[]

[Preconditioning]
 [./SMP]
   type = SMP
   full = false # true
  #  this will use a lot of memory. Could you try "full = false"? Or at least customize coupling elements.
 [../]
[]

[Executioner]
  type = Transient
  scheme = bdf2 # Type of time integration (2nd order backward euler), defaults to 1st order backward euler
  solve_type = PJFNK # Preconditioned JFNK (default)
  # petsc_options_iname = '-pc_type'
  # petsc_options_value = 'asm'

  petsc_options_iname = '-pc_type -sub_pc_type'
  petsc_options_value = 'asm ilu'

  l_tol = 1.0e-4 # Relative tolerance for linear solves
  l_max_its = 30 # Max number of linear iterations
  nl_max_its = 40 # 20 Max number of nonlinear iterations
  nl_rel_tol = 1.0e-10 # 1.0e-8 Absolute tolerance for nonlienar solves
  start_time = 0.0
  num_steps = 500
  dt = 0.0002

  [./TimeStepper]
    type = IterationAdaptiveDT
    cutback_factor = 0.9
    dt = 0.0002
    growth_factor = 1.1
    optimal_iterations = 8
  [../]

  # [./Adaptivity]
  #  initial_adaptivity = 4
  #  refine_fraction = 0.6
  #  coarsen_fraction = 0.1
  #  max_h_level = 4
  #  #  print_changed_info = true
  # [../]
[]

[Outputs]
  # file_base = sub1/sub2/subdir_output_out
  nemesis = true
  checkpoint = true
  csv = true
  [./console]
    type = Console
  [../]
[]
