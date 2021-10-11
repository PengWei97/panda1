# This simulation tests the TricrystalTripleJunctionIC

[Mesh]
  # Mesh block. Meshes can be read in or automatically generated.
  type = GeneratedMesh
  dim = 2 # Problem dimension
  nx = 110 # Number of elements in the x direction
  ny = 110 # Number of elements in the y direction
  xmax = 1000 # Maximum x-coordinate of mesh
  xmin = 0 # Minimum x-coordinate of mesh
  ymax = 1000 # Maximum y-coordinate of mesh
  ymin = 0 # Minimum y-coordinate of mesh
  elem_type = QUAD4 # Type of elements used in the mesh
  uniform_refine = 3
[]

[GlobalParams]
  # Parameters used by several kernels that are defined globally to simplify input file
  op_num = 3 # Number of order parameters used
  var_name_base = gr # base name of grains
  v = 'gr0 gr1 gr2' # Names of the grains
  theta1 = 135 # Angle the first grain makes at the triple junction
  theta2 = 100 # Angle the second grain makes at the triple junction
  length_scale = 1.0e-9 # Length scale in nm
  time_scale = 1.0e-9 # Time scale in ns
[]

[Variables]
  [./gr0]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
       type = TricrystalTripleJunctionIC
       op_index = 1
    [../]
  [../]

  [./gr1]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
       type = TricrystalTripleJunctionIC
       op_index = 2
    [../]
  [../]

  [./gr2]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
       type = TricrystalTripleJunctionIC
       op_index = 3
    [../]
  [../]
[]

[AuxVariables]
  [./bnds]
    # Variable used to visualize the grain boundaries in the simulation
    order = FIRST
    family = LAGRANGE
  [../]
  [./local_energy]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[Kernels]
  # Kernels block where the kernels defining the residual equations are set up
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
  [./local_free_energy]
    type = TotalFreeEnergy
    f_name = F
    variable = local_energy
    kappa_names = 'kappa_op kappa_op kappa_op'
    interfacial_vars = 'gr0 gr1 gr2'
  [../]
[]

[Materials]
  [./material]
    # Material properties
    type = GBEvolution
    T = 450 # Constant temperature of the simulation (for mobility calculation)
    wGB = 14 # Width of the diffuse GB
    GBmob0 = 2.5e-6 #m^4(Js) for copper from Schoenfelder1997
    Q = 0.23 #eV for copper from Schoenfelder1997
    GBenergy = 0.708 #J/m^2 from Schoenfelder1997
  [../]
  [./free_energy]
    type = DerivativeParsedMaterial
    f_name= F
    args = 'gr0 gr1 gr2'
    material_property_names = 'mu gamma_asymm'
    function = 'mu*( gr0^4/4.0 - gr0^2/2.0 + gr1^4/4.0 - gr1^2/2.0 + gr2^4/4.0 - gr2^2/2.0 + gamma_asymm*gr0^2*gr1^2  + gamma_asymm*gr0^2*gr2^2 + gamma_asymm*gr1^2*gr2^2) + 1.0/4.0'
    derivative_order = 2
    enable_jit = true
  [../]
[]

[Postprocessors]
  # Scalar postprocessors
  [./grain_tracker]
    type = FauxGrainTracker
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
  nl_max_its = 40
  nl_rel_tol = 1.0e-7
  start_time = 0.0
  num_steps = 50
  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 1.5
    growth_factor = 1.2
    cutback_factor = 0.8
    optimal_iterations = 8
  [../]
  [./Adaptivity]
    initial_adaptivity = 5
    refine_fraction = 0.8
    coarsen_fraction = 0.05
    max_h_level = 10
  [../]

[]

[Outputs]
  exodus = true # Outputs to the Exodus file format
  execute_on = 'initial timestep_end'
[]

# [Problem]
#   solve = false
# []
