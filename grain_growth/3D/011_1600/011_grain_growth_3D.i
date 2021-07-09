# This simulation predicts GB migration of a 2D copper polycrystal with 100 grains represented with 18 order parameters
# Mesh adaptivity and time step adaptivity are used
# An AuxVariable is used to calculate the grain boundary locations
# Postprocessors are used to record time step and the number of grains
# 001:开启自适应模式，可以正常演化，但是感觉演化不明显，w_GB = 125nm
# 002:基于001，添加输出晶粒体积，w_GB = 100nm,晶界正常演化，演化速率快于001
# 003:基于002，w_GB = 75
# 004：基于002和3D_6000_gr.i,length_scale = 1e-6,time_scale = 1,第一步不收敛
# 005：基于004，time_scale = 1e-3
# 006：基于004，w_GB = 100，xmax = 1000
# 007：基于006，w_GB= 3.5,xmax = 35，grain_raduis = 7.4255，time_scale = 1e-6
# 008: 基于007，扩大晶粒数目 grain_num = 200, xmax = 70

[Mesh]
  # Mesh block.  Meshes can be read in or automatically generated
  type = GeneratedMesh
  dim = 3 # Problem dimension
  nx = 280# Number of elements in the x-direction
  ny = 280 # Number of elements in the y-direction
  nz = 280
  xmin = 0    # minimum x-coordinate of the mesh
  xmax = 280 # maximum x-coordinate of the mesh
  ymin = 0    # minimum y-coordinate of the mesh
  ymax = 280 # maximum y-coordinate of the mesh
  zmin = 0
  zmax = 280
#   uniform_refine = 1 # Initial uniform refinement of the mesh

  parallel_type = distributed
[]

[GlobalParams]
  # Parameters used by several kernels that are defined globally to simplify input file
  op_num = 30 # Number of order parameters used
  var_name_base = gr # Base name of grains
  order = CONSTANT
  family = MONOMIAL
[]

[Variables]
  # Variable block, where all variables in the simulation are declared
  [./PolycrystalVariables]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[UserObjects]
  [./voronoi]
    type = PolycrystalVoronoi
    grain_num = 12800 # Number of grains
    rand_seed = 6000
    coloring_algorithm = jp
  [../]
  [./grain_tracker]
    type = GrainTracker
    compute_var_to_feature_map = true
    threshold = 0.2
    connecting_threshold = 0.08
    compute_halo_maps = true # Only necessary for displaying HALOS
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
  # Dependent variables
  [./bnds]
    # Variable used to visualize the grain boundaries in the simulation
    order = FIRST
    family = LAGRANGE
  [../]
  [./unique_grains]
  [../]

  # [./var_indices]
  # [../]

  # [./ghost_regions]
  # [../]

  # [./halos]
  # [../]

  # [./halo0]
  # [../]

  # [./halo1]
  # [../]

  # [./halo2]
  # [../]

  # [./halo3]
  # [../]

  # [./halo4]
  # [../]

  # [./halo5]
  # [../]

  # [./halo6]
  # [../]

  # [./halo7]
  # [../]

  # [./halo8]
  # [../]

  # [./halo9]
  # [../]

  # [./halo10]
  # [../]

  # [./halo11]
  # [../]

  # [./halo12]
  # [../]

  # [./halo13]
  # [../]

  # [./halo14]
  # [../]

  # [./proc]
  # [../]
[]

[Kernels]
  # Kernel block, where the kernels defining the residual equations are set up.
  [./PolycrystalKernel]
    # Custom action creating all necessary kernels for grain growth.  All input parameters are up in GlobalParams
  [../]
[]

[AuxKernels]
  # AuxKernel block, defining the equations used to calculate the auxvars
  [./bnds_aux]
    # AuxKernel that calculates the GB term
    type = BndsCalcAux
    variable = bnds
    execute_on = 'initial timestep_end'
  [../]
  [./unique_grains]
    type = FeatureFloodCountAux
    variable = unique_grains
    flood_counter = grain_tracker
    field_display = UNIQUE_REGION
    execute_on = 'initial timestep_end'
  [../]
  # [./var_indices]
  #   type = FeatureFloodCountAux
  #   variable = var_indices
  #   flood_counter = grain_tracker
  #   field_display = VARIABLE_COLORING
  #   execute_on = 'initial timestep_end'
  # [../]
  # [./ghosted_entities]
  #   type = FeatureFloodCountAux
  #   variable = ghost_regions
  #   flood_counter = grain_tracker
  #   field_display = GHOSTED_ENTITIES
  #   execute_on = 'initial timestep_end'
  # [../]
  # [./halos]
  #   type = FeatureFloodCountAux
  #   variable = halos
  #   flood_counter = voronoi
  #   field_display = HALOS
  #   execute_on = 'initial timestep_end'
  # [../]
  # [./halo0]
  #   type = FeatureFloodCountAux
  #   variable = halo0
  #   map_index = 0
  #   field_display = HALOS
  #   flood_counter = grain_tracker
  #   execute_on = 'initial timestep_end'
  # [../]
  # [./halo1]
  #   type = FeatureFloodCountAux
  #   variable = halo1
  #   map_index = 1
  #   field_display = HALOS
  #   flood_counter = grain_tracker
  # [../]
  # [./halo2]
  #   type = FeatureFloodCountAux
  #   variable = halo2
  #   map_index = 2
  #   field_display = HALOS
  #   flood_counter = grain_tracker
  # [../]
  # [./halo3]
  #   type = FeatureFloodCountAux
  #   variable = halo3
  #   map_index = 3
  #   field_display = HALOS
  #   flood_counter = grain_tracker
  # [../]
  # [./halo4]
  #   type = FeatureFloodCountAux
  #   variable = halo4
  #   map_index = 4
  #   field_display = HALOS
  #   flood_counter = grain_tracker
  # [../]
  # [./halo5]
  #   type = FeatureFloodCountAux
  #   variable = halo5
  #   map_index = 5
  #   field_display = HALOS
  #   flood_counter = grain_tracker
  # [../]
  # [./halo6]
  #   type = FeatureFloodCountAux
  #   variable = halo6
  #   map_index = 6
  #   field_display = HALOS
  #   flood_counter = grain_tracker
  # [../]
  # [./halo7]
  #   type = FeatureFloodCountAux
  #   variable = halo7
  #   map_index = 7
  #   field_display = HALOS
  #   flood_counter = grain_tracker
  # [../]
  # [./halo8]
  #   type = FeatureFloodCountAux
  #   variable = halo8
  #   map_index = 8
  #   field_display = HALOS
  #   flood_counter = grain_tracker
  # [../]
  # [./halo9]
  #   type = FeatureFloodCountAux
  #   variable = halo9
  #   map_index = 9
  #   field_display = HALOS
  #   flood_counter = grain_tracker
  # [../]
  # [./halo10]
  #   type = FeatureFloodCountAux
  #   variable = halo10
  #   map_index = 10
  #   field_display = HALOS
  #   flood_counter = grain_tracker
  # [../]
  # [./halo11]
  #   type = FeatureFloodCountAux
  #   variable = halo11
  #   map_index = 11
  #   field_display = HALOS
  #   flood_counter = grain_tracker
  # [../]
  # [./halo12]
  #   type = FeatureFloodCountAux
  #   variable = halo12
  #   map_index = 12
  #   field_display = HALOS
  #   flood_counter = grain_tracker
  # [../]
  # [./halo13]
  #   type = FeatureFloodCountAux
  #   variable = halo13
  #   map_index = 13
  #   field_display = HALOS
  #   flood_counter = grain_tracker
  # [../]
  # [./halo14]
  #   type = FeatureFloodCountAux
  #   variable = halo14
  #   map_index = 14
  #   field_display = HALOS
  #   flood_counter = grain_tracker
  # [../]
  # [./proc]
  #   type = ProcessorIDAux
  #   variable = proc
  #   execute_on = 'initial timestep_end'
  # [../]

[]

[Materials]
  [./CuGrGr]
    # Material properties
    type = GBEvolution
    T = 450 # Constant temperature of the simulation (for mobility calculation)
    wGB = 3.5 # Width of the diffuse GB μm
    GBmob0 = 2.5e-6 #m^4(Js) for copper from Schoenfelder1997
    Q = 0.23 #eV for copper from Schoenfelder1997
    GBenergy = 0.708 #J/m^2 from Schoenfelder1997

    # wGB = 3 # um
    length_scale = 1.0e-6 # μm
    time_scale = 1.0e-6 # μs
  [../]
[]

[Postprocessors]
  # Scalar postprocessors
  [./ngrains]
    type = FeatureFloodCount
    variable = bnds
    threshold = 0.7
  [../]
  [./dt]
    # Outputs the current time step
    type = TimestepSize
  [../]
  [./run_time]
    type = PerfGraphData
    section_name = "Root"
    data_type = total
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
  type = Transient # Type of executioner, here it is transient with an adaptive time step
  scheme = bdf2 # Type of time integration (2nd order backward euler), defaults to 1st order backward euler

  #Preconditioned JFNK (default)
  solve_type = 'PJFNK'

  # Uses newton iteration to solve the problem.
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'asm'

  l_max_its = 30 # Max number of linear iterations
  l_tol = 1e-4 # Relative tolerance for linear solves
  nl_max_its = 40 # Max number of nonlinear iterations
  nl_rel_tol = 1e-10 # Absolute tolerance for nonlienar solves

  start_time = 0.0
  end_time = 80000

  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 25 # Initial time step.  In this simulation it changes.
    optimal_iterations = 6 # Time step will adapt to maintain this number of nonlinear iterations
  [../]

#  [./Adaptivity]
#    # Block that turns on mesh adaptivity. Note that mesh will never coarsen beyond initial mesh (before uniform refinement)
#    initial_adaptivity = 2 # Number of times mesh is adapted to initial condition
#    refine_fraction = 0.6 # Fraction of high error that will be refined
#    coarsen_fraction = 0.1 # Fraction of low error that will coarsened
#    max_h_level = 3 # Max number of refinements used, starting from initial mesh (before uniform refinement)
#  [../]
[]

[Outputs]
  [./csv]
    type = CSV
    interval = 1
  [../]
  [./exodus1]
    file_base = 3D_200_0-10000
    type = Exodus
    interval = 1
#     append_date = true # 添加输出时间
    start_time = 0
    end_time = 10000
#     end_step = 5 # 在计算步为第五步的时候停止
  [../]
  [./exodus2]
    file_base = 3D_200_10000-20000
    type = Exodus
    interval = 1
    start_time = 10000
    end_time = 20000
  [../]
  [./exodus3]
    file_base = 3D_200_20000-30000
    type = Exodus
    interval = 1
    append_date = true # 添加输出时间
    start_time = 20000
    end_time = 30000
  [../]
  [./exodus4]
    file_base = 3D_200_30000-40000
    type = Exodus
    interval = 1
    append_date = true # 添加输出时间
    start_time = 30000
    end_time = 40000
  [../]
    [./exodus5]
    file_base = 3D_200_40000-50000
    type = Exodus
    interval = 1
    append_date = true # 添加输出时间
    start_time = 40000
    end_time = 5000
  [../]
  [./exodus6]
    file_base = 3D_200_50000-60000
    type = Exodus
    interval = 1
    append_date = true # 添加输出时间
    start_time = 50000
    end_time = 60000
  [../]
  [./exodus7]
    file_base = 3D_200_60000-70000
    type = Exodus
    interval = 1
    append_date = true # 添加输出时间
    start_time = 60000
    end_time = 70000
  [../]
  [./exodus8]
    file_base = 3D_200_70000-80000
    type = Exodus
    interval = 1
    append_date = true # 添加输出时间
    start_time = 70000
    end_time = 80000
  [../]
  [./perf_graph]
    type = PerfGraphOutput
    execute_on = 'initial final'  # Default is "final"
    level = 7                     # Default is 1
  [../]
  [out]
    type = Checkpoint
    interval = 10
    num_files = 6
  []
[]
