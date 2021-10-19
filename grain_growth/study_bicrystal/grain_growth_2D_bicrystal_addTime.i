# This simulation predicts GB migration of a 2D copper polycrystal with 15 grains
# Mesh adaptivity and time step adaptivity are used
# An AuxVariable is used to calculate the grain boundary locations
# Postprocessors are used to record time step and the number of grains
# We are not using the GrainTracker in this example so the number
# of order paramaters must match the number of grains.
# the grain raduis of copper 200 μm, initial 50 μm
my_GBmob0 = 2.5e-6
my_length_scale = 1.0e-6
my_time_scale = 1.0e-3
# my_wGB = 4 #8 # μm
my_filename = 'addTime_4wGB'

[Mesh]
  # Mesh block.  Meshes can be read in or automatically generated
  type = GeneratedMesh
  dim = 2
  nx = 50
  ny = 50
  nz = 0
  xmax = 500
  ymax = 500
  elem_type = QUAD4

  uniform_refine = 2
[]

[GlobalParams]
  # Parameters used by several kernels that are defined globally to simplify input file
  op_num = 2 # Number of grains
  var_name_base = gr # Base name of grains
[]

# [UserObjects]
#   [./voronoi]
#     type = PolycrystalVoronoi
#     grain_num = 2
#     rand_seed = 42
#     coloring_algorithm = jp # We must use bt to force the UserObject to assign one grain to each op
#   [../]
# []

[ICs]
  [./PolycrystalICs]
    [./BicrystalCircleGrainIC]
      radius = 200
      x = 250
      y = 250
      int_width = 10
    [../]
  [../]
[]

[Variables]
  # Variable block, where all variables in the simulation are declared
  [./PolycrystalVariables]
    # Custom action that created all of the grain variables
    order = FIRST # element type used by each grain variable
    family = LAGRANGE
  [../]
[]

[AuxVariables]
#active = ''
  # Dependent variables
  [./bnds]
    # Variable used to visualize the grain boundaries in the simulation
    order = FIRST
    family = LAGRANGE
  [../]
  # [./unique_grains]
  #   order = CONSTANT
  #   family = MONOMIAL
  # [../]
  # [./var_indices]
  #   order = CONSTANT
  #   family = MONOMIAL
  # [../]
  [./local_energy]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[Kernels]
  # Kernel block, where the kernels defining the residual equations are set up.
  [./PolycrystalKernel]
    # Custom action creating all necessary kernels for grain growth.  All input parameters are up in GlobalParams
  [../]
[]

[AuxKernels]
#active = ''
  # AuxKernel block, defining the equations used to calculate the auxvars
  [./bnds_aux]
    # AuxKernel that calculates the GB term
    type = BndsCalcAux
    variable = bnds
    execute_on = timestep_end
  [../]
  # [./unique_grains]
  #   type = FeatureFloodCountAux
  #   variable = unique_grains
  #   execute_on = timestep_end
  #   flood_counter = grain_tracker
  #   field_display = UNIQUE_REGION
  # [../]
  # [./var_indices]
  #   type = FeatureFloodCountAux
  #   variable = var_indices
  #   execute_on = timestep_end
  #   flood_counter = grain_tracker
  #   field_display = VARIABLE_COLORING
  # [../]
  [./local_free_energy]
    type = TotalFreeEnergy
    f_name = F
    variable = local_energy
    kappa_names = 'kappa_op kappa_op'
    interfacial_vars = 'gr0 gr1'
  [../]
[]

[BCs]
  # Boundary Condition block
  [./Periodic]
    [./top_bottom]
      auto_direction = 'x y' # Makes problem periodic in the x and y directions
    [../]
  [../]
[]

[Materials]
  [./CuGrGr]
    # Material properties
    type = GBEvolution # Quantitative material properties for copper grain growth.  Dimensions are nm and ns
    GBmob0 = ${my_GBmob0} # 2.5e-6 #Mobility prefactor for Cu from Schonfelder1997
    GBenergy = 0.708 #GB energy for Cu from Schonfelder1997
    Q = 0.23 #Activation energy for grain growth from Schonfelder 1997
    T = 500 # K   #Constant temperature of the simulation (for mobility calculation)
    # wGB = ${my_wGB} # nm      #Width of the diffuse GB
    wGB = 8 # nm      #Width of the diffuse GB
    length_scale = ${my_length_scale}
    time_scale = ${my_time_scale}
    outputs = exodus
  [../]
  [./free_energy]
    type = DerivativeParsedMaterial
    f_name= F
    args = 'gr0 gr1'
    material_property_names = 'mu gamma_asymm'
    function = 'mu*( gr0^4/4.0 - gr0^2/2.0 + gr1^4/4.0 - gr1^2/2.0 + gamma_asymm*gr0^2*gr1^2 + 1.0/4.0) '
    derivative_order = 2
    enable_jit = true
  [../]
[]

[Postprocessors]
  # active = 'dt '
  # # Scalar postprocessors
  [./dt]
    # Outputs the current time step
    type = TimestepSize
  [../]
  [./dofs]
    type = NumDOFs
  [../]
  [./run_time]
    type = PerfGraphData
    section_name = "Root"
    data_type = total
  [../]
  [./gr0_area]
    type = ElementIntegralVariablePostprocessor
    variable = gr0
  [../]
[]

[Executioner]
  type = Transient # Type of executioner, here it is transient with an adaptive time step
  scheme = bdf2 # Type of time integration (2nd order backward euler), defaults to 1st order backward euler

  #Preconditioned JFNK (default)
  solve_type = 'PJFNK'

  petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart -mat_mffd_type'
  petsc_options_value = 'hypre    boomeramg      101                ds'

  l_max_its = 30 # Max number of linear iterations
  l_tol = 1e-4 # Relative tolerance for linear solves
  nl_max_its = 40 # Max number of nonlinear iterations
  nl_abs_tol = 1e-11 # Relative tolerance for nonlienar solves
  nl_rel_tol = 1e-8 # Absolute tolerance for nonlienar solves

  start_time = 0.0
  end_time = 4000

  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 1.5
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
#   file_base = ./grain_growth_2D_bicrystal/bicrystal
  file_base = ./${my_filename}/out_${my_filename}
  exodus = true
  csv = true
  [./my_checkpoint]
    type = Checkpoint
    num_files = 6
    interval = 5
  [../]
  [./console]
    type = Console
    max_rows = 20
  [../]
  [pgraph]
    type = PerfGraphOutput
    execute_on = 'initial final'  # Default is "final"
    level = 2                     # Default is 1
    heaviest_branch = true        # Default is false
    heaviest_sections = 7         # Default is 0
  []
[]
