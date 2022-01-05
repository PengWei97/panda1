my_filename = 'Pattern1_gg_polycrystal' 
my_num_adaptivity = 3
my_interval = 2

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 25
  ny = 25
  xmax = 1000
  ymax = 1000
  elem_type = QUAD4
  # uniform_refine = 2
[]

[GlobalParams]
  # CahnHilliard needs the third derivatives
  op_num = 2
  var_name_base = gr
  derivative_order = 2
  # enable_jit = true # 
  displacements = 'disp_x disp_y'
[]


[ICs]
  [./PolycrystalICs]
    [./BicrystalCircleGrainIC]
      radius = 300
      x = 500
      y = 500
      int_width = 10
    [../]
  [../]
[]

# AuxVars to compute the free energy density for outputting
[AuxVariables]
  # [./local_energy]
  #   order = CONSTANT
  #   family = MONOMIAL
  # [../]
  # [./cross_energy]
  #   order = CONSTANT
  #   family = MONOMIAL
  # [../]
  [./bnds]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[AuxKernels]
  # [./local_free_energy]
  #   type = TotalFreeEnergy
  #   variable = local_energy
  #   interfacial_vars = 'c'
  #   kappa_names = 'kappa_c'
  #   additional_free_energy = cross_energy
  # [../]
  # [./cross_terms]
  #   type = CrossTermGradientFreeEnergy
  #   variable = cross_energy
  #   interfacial_vars = 'gr0 gr1 eta3'
  #   kappa_names = 'kappa11 kappa12 kappa13
  #                  kappa21 kappa22 kappa23
  #                  kappa31 kappa32 kappa33'
  # [../]
  [./BndsCalc]
    type = BndsCalcAux
    variable = bnds
  [../]
[]

[Variables]
  # Order parameter for the Matrix
  [./PolycrystalVariables]
    # gr0 gr1
  [../]

  # Mesh displacement
  [./disp_x]
    order = FIRST
    family = LAGRANGE
  [../]
  [./disp_y]
    order = FIRST
    family = LAGRANGE
  [../]

  # Lagrange-multiplier
  # [./lambda] 
  #   order = FIRST
  #   family = LAGRANGE
  #   initial_condition = 1.0
  # [../]
[]

[UserObjects]
  [./grain_tracker]
    type = GrainTracker
    threshold = 0.2
    connecting_threshold = 0.08
    compute_halo_maps = true # Only necessary for displaying HALOS
  [../]
[]

[Kernels]
  # Set up stress divergence kernels
  [./TensorMechanics]
  [../]
  # [./PolycrystalElasticDrivingForce]
  # [../]
  [./PolycrystalKernel]
  [../]
  [./ACBulk1]
    type = AllenCahn
    variable = gr0
    args = 'gr1'
    mob_name = L
    f_name = F
  [../]
  # [./lagrange1]
  #   type = SwitchingFunctionConstraintEta # Lagrange multiplier kernel to constrain the sum of all switching functions in a multiphase system.
  #   # the total weight of all phase free energy contributions at each point in the simulation volume is exactly unity
  #   variable = gr0
  #   h_name   = h1 # Switching Function Materials that provides h(eta_i)
  #   lambda = lambda # Lagrange multiplier
  # [../]

  # Allen-Cahn and Lagrange-multiplier constraint kernels for order parameter 2
  [./ACBulk2]
    type = AllenCahn
    variable = gr1
    args = 'gr0'
    mob_name = L
    f_name = F
  [../]
  # [./lagrange2]
  #   type = SwitchingFunctionConstraintEta
  #   variable = gr1
  #   h_name   = h2
  #   lambda = lambda
  # [../]

  # Lagrange-multiplier constraint kernel for lambda
  # [./lagrange]
  #   type = SwitchingFunctionConstraintLagrange
  #   variable = lambda
  #   etas    = 'gr0 gr1'
  #   h_names = 'h1   h2'
  #   epsilon = 1e-6
  # [../]
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

  # We use this to output the level of constraint enforcement
  # ideally it should be 0 everywhere, if the constraint is fully enforced
  # [./etasummat]
  #   type = ParsedMaterial
  #   f_name = etasum
  #   args = 'gr0 gr1'
  #   material_property_names = 'h1 h2'
  #   function = 'h1+h2-1'
  #   outputs = my_exodus
  # [../]

  # matrix phase
  [./elasticity_tensor_1]
    type = ComputeElasticityTensor
    base_name = phase1
    C_ijkl = '1.27e5 0.708e5 0.708e5 1.27e5 0.708e5 1.27e5 0.7355e5 0.7355e5 0.7355e5' # FCC Cu
    fill_method = symmetric9
    euler_angle_1 = 0. # same as above but opposite direction because _transpose_ gets built from these angles
    euler_angle_2 = 0.
    euler_angle_3 = 0.
  [../]
  [./strain_1]
    type = ComputeSmallStrain
    base_name = phase1
    displacements = 'disp_x disp_y'
  [../]
  [./stress_1]
    type = ComputeLinearElasticStress
    base_name = phase1
  [../]

  # oversized phase
  [./elasticity_tensor_2]
    type = ComputeElasticityTensor
    base_name = phase2
    C_ijkl = '1.27e5 0.708e5 0.708e5 1.27e5 0.708e5 1.27e5 0.7355e5 0.7355e5 0.7355e5' # FCC Cu
    fill_method = symmetric9
    euler_angle_1 = 45. # same as above but opposite direction because _transpose_ gets built from these angles
    euler_angle_2 = 0.
    euler_angle_3 = 0.
  [../]
  [./strain_2]
    type = ComputeSmallStrain
    base_name = phase2
    displacements = 'disp_x disp_y'
    # eigenstrain_names = eigenstrain
  [../]
  [./stress_2]
    type = ComputeLinearElasticStress
    base_name = phase2
  [../]

  # switching functions
  [./switching1]
    type = SwitchingFunctionMaterial
    function_name = h1
    eta = gr0
    h_order = SIMPLE
  [../]
  [./switching2]
    type = SwitchingFunctionMaterial
    function_name = h2
    eta = gr1
    h_order = SIMPLE
  [../]

  # elastic free energies
  [./elastic_free_energy_1]
    type = ElasticEnergyMaterial
    base_name = phase1
    f_name = Fe1
    derivative_order = 2
    args = 'gr0 gr1' # should be empty
  [../]
  [./elastic_free_energy_2]
    type = ElasticEnergyMaterial
    base_name = phase2
    f_name = Fe2
    derivative_order = 2
    args = 'gr0 gr1' # should be empty
  [../]

  # # # phase free energies (chemical + elastic)
  # [./phase_free_energy_1]
  #   type = DerivativeSumMaterial
  #   f_name = F1
  #   sum_materials = 'Fc1 Fe1'
  #   args = 'gr0 gr1'
  #   derivative_order = 2
  # [../]
  # [./phase_free_energy_2]
  #   type = DerivativeSumMaterial
  #   f_name = F2
  #   sum_materials = 'Fc2 Fe2'
  #   args = 'gr0 gr1'
  #   derivative_order = 2
  # [../]

  [./barrier]
    type = MultiBarrierFunctionMaterial
    etas = 'gr0 gr1'
    # Double well phase transformation barrier free energy contribution.
  [../]
  
  # global free energy
  [./free_energy]
    type = DerivativeMultiPhaseMaterial
    # Two phase material that combines n phase materials using a switching function with and n non-conserved order parameters
    f_name = F
    fi_names = 'Fe1  Fe2'
    hi_names = 'h1  h2'
    etas     = 'gr0 gr1'
    # args = 'gr0 gr1' # Arguments of the fi free energies - use vector coupling
    # 如果不求c可以不输入args
    W = 0 # Energy barrier for the phase transformation from A to B
  [../]

  # Generate the global stress from the phase stresses
  [./global_stress]
    type = MultiPhaseStressMaterial
    phase_base = 'phase1 phase2'
    h          = 'h1     h2'
  [../]
[]

[BCs]
  # the boundary conditions on the displacement enforce periodicity
  # at zero total shear and constant volume
  [./Periodic]
    [./All]
      auto_direction = 'x'
      variable = 'gr0 gr1'
    [../]
  [../]
  [./bottom_y]
    type = DirichletBC
    variable = disp_y
    boundary = 'bottom'
    value = 0
  [../]
  [./top_y]
    type = DirichletBC
    variable = disp_y
    boundary = 'top'
    value = 0
  [../]
  [./left_x]
    type = DirichletBC
    variable = disp_x
    boundary = 'left'
    value = 0
  [../]
  [./right_x]
    type = DirichletBC
    variable = disp_x
    boundary = 'right'
    value = 0
  [../]

  # [./Periodic]
  #   [./disp_x]
  #     auto_direction = 'y'
  #   [../]
  #   [./disp_y]
  #     auto_direction = 'x'
  #   [../]

  #   # all other phase field variables are fully periodic

  #   [./lambda]
  #     auto_direction = 'x y'
  #   [../]
  # [../]
[]

[Preconditioning]
  [./SMP]
    type = SMP
    full = true
  [../]
[]

# We monitor the total free energy and the total solute concentration (should be constant)
[Postprocessors]
  # [./total_free_energy]
  #   type = ElementIntegralVariablePostprocessor
  #   variable = local_energy
  # [../]
  # [./total_solute]
  #   type = ElementIntegralVariablePostprocessor
  #   variable = c
  # [../]
  [./gr1area]
    type = ElementIntegralVariablePostprocessor
    variable = gr0
  [../]
[]

[Executioner]
  type = Transient
  scheme = bdf2

  solve_type = 'PJFNK'
  petsc_options_iname = '-pc_type -sub_pc_type'
  petsc_options_value = 'asm      ilu'

  l_max_its = 30
  nl_max_its = 10
  l_tol = 1.0e-4
  nl_rel_tol = 1.0e-8
  nl_abs_tol = 1.0e-10
  start_time = 0.0
  num_steps = 200

  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 1.5
    growth_factor = 1.2
    cutback_factor = 0.8
    optimal_iterations = 8
  [../]
  [./Adaptivity]
    initial_adaptivity = ${my_num_adaptivity}
    refine_fraction = 0.8 # The fraction of elements or error to refine. Should be between 0 and 1.
    coarsen_fraction = 0.05 # 0.05 Fraction of low error that will coarsened
    max_h_level = ${my_num_adaptivity}
  [../]
[]

[Outputs]
  file_base = ./${my_filename}/out_${my_filename}
  # execute_on = 'timestep_end'
  [./my_exodus]
    type = Exodus
    interval = ${my_interval} # The interval at which time steps are output
    # sync_times = '10 50 100 500 1000 5000 10000 50000 100000'
    # sync_only = true
    sequence = true
  [../]
  
  [./csv]
    type = CSV
    # interval = 1
  [../]
  [pgraph]
    type = PerfGraphOutput
    execute_on = 'initial final'  # Default is "final"
    level = 2                     # Default is 1
    heaviest_branch = true        # Default is false
    heaviest_sections = 7         # Default is 0
  []  
[]

[Debug]
  # show_var_residual_norms = true
[]
