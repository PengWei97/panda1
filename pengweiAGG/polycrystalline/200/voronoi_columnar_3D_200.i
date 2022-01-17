my_filename = 'GG_exmaples200_3D_statistic' 

[Mesh]
  type = GeneratedMesh
  dim = 3
  nx = 40
  ny = 40
  nz = 1
  xmin = 0
  xmax = 2000
  ymin = 0
  ymax = 2000
  zmin = 0
  zmax = 25
  elem_type = HEX8
  uniform_refine = 1

  parallel_type = distributed # Periodic BCs
[]

[GlobalParams]
  op_num = 15
  var_name_base = gr
  grain_num = 200
[]

[Variables]
  [./PolycrystalVariables]
  [../]
[]

[UserObjects]
  [./voronoi]
    type = PolycrystalVoronoi
    rand_seed = 300
    columnar_3D = true
    coloring_algorithm = jp
  [../]
  [./grain_tracker]
    type = GrainTracker
    threshold = 0.2
    connecting_threshold = 0.08
    compute_halo_maps = true # Only necessary for displaying HALOS
    compute_var_to_feature_map = true  # for FeatureVolumeVectorPostprocessor
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
[]

[Kernels]
  [./PolycrystalKernel]
  [../]
[]

[AuxKernels]
  [./BndsCalc]
    type = BndsCalcAux
    variable = bnds
    execute_on = timestep_end
  [../]
  [./unique_grains]
    type = FeatureFloodCountAux
    variable = unique_grains
    flood_counter = grain_tracker
    field_display = UNIQUE_REGION
    execute_on = 'initial timestep_end'
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
    # type = GBEvolution
    # T = 500 # K
    # wGB = 60 # nm
    # GBmob0 = 2.5e-6 #m^4/(Js) from Schoenfelder 1997
    # Q = 0.23 #Migration energy in eV
    # GBenergy = 0.708 #GB energy in J/m^2

    type = GBEvolution
    block = 0
    T = 500 # K
    wGB = 15 # nm
    GBmob0 = 2.5e-6 # m^4/(Js) from Schoenfelder 1997
    Q = 0.23 # Migration energy in eV
    GBenergy = 0.708 # GB energy in J/m^2
    # length_scale = 1.0e-6
    # time_scale = 1.0
  [../]
[]

[Postprocessors]
  # active = ''
  [./ngrains]
    type = FeatureFloodCount
    variable = bnds
    threshold = 0.7
  [../]
  [./dofs]
    type = NumDOFs
  [../]
  [./run_time]
    type = PerfGraphData
    section_name = "Root"
    data_type = total
  [../]
  [./dt]
    type = TimestepSize
  [../]
  [./avg_grain_volumes]
    type = AverageGrainVolume
    feature_counter = grain_tracker
    execute_on = 'initial timestep_end'
  [../]
[]

[VectorPostprocessors]
  [./grain_volumes]
    type = FeatureVolumeVectorPostprocessor
    flood_counter = grain_tracker
    execute_on = 'initial timestep_end'
  #  output_centroids = true
  [../]
[]

[Executioner]
  type = Transient
  scheme = bdf2
  solve_type = PJFNK
  petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart -pc_hypre_boomeramg_strong_threshold'
  petsc_options_value = 'hypre boomeramg 31 0.7'

  l_tol = 1.0e-4
  l_max_its = 30
  nl_max_its = 25
  nl_rel_tol = 1.0e-7

  start_time = 0.0
  end_time = 10000.0

  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 1.5
    growth_factor = 1.2
    cutback_factor = 0.8
    optimal_iterations = 8
  [../]

  [./Adaptivity]
    initial_adaptivity = 4
    refine_fraction = 0.8
    coarsen_fraction = 0.1
    max_h_level = 4
  [../]
[]

[Outputs]
  file_base = ./${my_filename}/out_${my_filename}
  # exodus = true
  [./csv]
    type = CSV
    # interval = 1
  [../]
  [./exodus]
    type = Nemesis
    interval = 2
  [../]
  [pgraph]
    type = PerfGraphOutput
    execute_on = 'initial timestep_end'  # Default is "final"
    level = 2                     # Default is 1
    heaviest_branch = true        # Default is false
    heaviest_sections = 7         # Default is 0
  []
[]
