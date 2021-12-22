
my_filename = 'test01_elastic_multiphase_02'
my_function = 'if(t<4,0.1*t,0.4+0.002*sin(10*pi*t))' # 0.001s^{-1}
my_end_time = 60 #60 #4e2

my_xmax = 30 # 3.0e3 30.0 #
my_ymax = 10 # 1.0e3 10.0 #

my_radius = 10.0 #10.0e1 # 1.0e3
my_nx = 100 # 200 # 200 # 50
my_ny = 25 # 50 # 50 # 20

my_time_scale = 1.0e-9
my_length_scale = 1.0e-9

my_wGB = 1.5 # 8 for 002

[Functions]
  [./timestep_fn]
    type = ParsedFunction
    value = 'if(t<4,0.1,5)'
  [../]
[]

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = ${my_nx}
  ny = ${my_ny}
  xmax = ${my_xmax}
  ymax = ${my_ymax}
  elem_type = QUAD4
  # uniform_refine = 2
[]

[GlobalParams]
  op_num = 2
  var_name_base = gr
  displacements = 'disp_x disp_y'
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
  [./stress11]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress12]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress22]
    order = CONSTANT
    family = MONOMIAL
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
  #   interfacial_vars = 'eta1 eta2 eta3'
  #   kappa_names = 'kappa11 kappa12 kappa13
  #                  kappa21 kappa22 kappa23
  #                  kappa31 kappa32 kappa33'
  # [../]
  [./stress11]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress11
    index_i = 0
    index_j = 0
  [../]
  [./stress12]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress12
    index_i = 0
    index_j = 1
  [../]
  [./stress22]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress22
    index_i = 1
    index_j = 1
  [../]
[]

[Variables]

  # Order parameter for the Matrix
  [./gr0]
    order = FIRST
    family = LAGRANGE
  [../]
  # Order parameters for the 2 different inclusion orientations
  [./gr1]
    order = FIRST
    family = LAGRANGE
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
[]

[ICs]
  [./PolycrystalICs]
    [./BicrystalBoundingBoxIC]
      x1 = 0
      y1 = 0
      x2 = ${my_radius}
      y2 = ${my_ymax}   
    [../]
  [../]
[]

[Kernels]
  # Set up stress divergence kernels
  [./TensorMechanics]
  [../]
  [./PolycrystalKernel]
  [../]
  # [./AC_ElasticDrivingForce_gr0]
  #   type = ACGrGrElasticDrivingForce
  #   D_tensor_name = delasticity_tensor/dgr0
  #   variable = gr0
  #   # D_stress = _D_elastic_tensor[_qp] * strain;
  #   # ACBulk.h
  # [../]
  # [./AC_ElasticDrivingForce_gr1]
  #   type = ACGrGrElasticDrivingForce
  #   variable = gr1
  #   D_tensor_name = delasticity_tensor/dgr1
  #   # D_stress = _D_elastic_tensor[_qp] * strain;
  # [../]
  [./ACBulk_ElasticEnergy_gr0]
    type = AllenCahn
    variable = gr0
    args = 'gr1'
    mob_name = L
    f_name = Fe1
  [../]
  [./ACBulk_ElasticEnergy_gr1]
    type = AllenCahn
    variable = gr1
    args = 'gr0'
    mob_name = L
    f_name = Fe2
  [../]
[]

# [Modules/TensorMechanics/Master]
#   [./all]
#     # strain = FINITE
#     # displacements = 'disp_x disp_y'  
#     # use_displaced_mesh = true
#     # add_variables = true
#     # use_automatic_differentiation = true
#     strain = SMALL # SMALL # FINITE
#     incremental = true
#   [../]
# []

[Materials]
  [./Copper]
    type = GBEvolution
    block = 0
    T = 500 # K
    wGB = ${my_wGB} #15 # nm
    GBmob0 = 2.5e-6 #m^4/(Js) from Schoenfelder 1997
    Q = 0.23 # Migration energy in eV
    GBenergy = 0.708 # GB energy in J/m^2
    time_scale = ${my_time_scale}
    length_scale = ${my_length_scale}
    # 提供材料参数
  [../]

  # gr0
  [./elasticity_tensor_1]
    type = ComputeElasticityTensor
    base_name = phase1
    C_ijkl = '1.684e5 1.214e5 1.214e5 1.684e5 1.214e5 1.684e5 0.75e5 0.75e5 0.75e5'
    fill_method = symmetric9
    euler_angle_1 = 0.0
    euler_angle_2 = 0.0
    euler_angle_3 = 0.0
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
    C_ijkl = '1.684e5 1.214e5 1.214e5 1.684e5 1.214e5 1.684e5 0.75e5 0.75e5 0.75e5'
    fill_method = symmetric9
    euler_angle_1 = 45.0
    euler_angle_2 = 0.0
    euler_angle_3 = 0.0
  [../]
  [./strain_2]
    type = ComputeSmallStrain
    base_name = phase2
    displacements = 'disp_x disp_y'
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


  # #   # chemical free energies
  # [./chemical_free_energy_1]
  #   type = DerivativeParsedMaterial
  #   f_name = Fc1
  #   material_property_names = 'mu gamma_asymm'
  #   function = 'mu*( gr0^4/4.0 - gr0^2/2.0 + gr1^4/4.0 - gr1^2/2.0 + gamma_asymm*gr0^2*gr1^2) + 1.0/4.0'
  #   args = 'gr0 gr1' # use vector coupling
  #   derivative_order = 2
  #   enable_jit = true
  # [../]
  # [./chemical_free_energy_2]
  #   type = DerivativeParsedMaterial
  #   f_name = Fc2
  #   function = '1/4*gr1^4-1/2*gr1^2+gamma_asymm*gr0^2*gr1^2'
  #   args = 'gr1' # use vector coupling
  #   derivative_order = 2
  # [../]

  # [./barrier]
  #   type = MultiBarrierFunctionMaterial
  #   etas = 'eta1 eta2 eta3'
  # [../]

  # elastic free energies
  [./elastic_free_energy_1]
    type = ElasticEnergyMaterial
    base_name = phase1
    f_name = Fe1
    derivative_order = 2
    args = 'gr0' # should be empty
  [../]
  [./elastic_free_energy_2]
    type = ElasticEnergyMaterial
    base_name = phase2
    f_name = Fe2
    derivative_order = 2
    args = 'gr1' # should be empty
  [../]

  # global free energy
  [./elastic_free_energy]
    type = DerivativeMultiPhaseMaterial
    f_name = F # base name of the free energy function (used to name the material properties)
    fi_names = 'Fe1  Fe1' # List of free energies for the n phases
    hi_names = 'h1  h2' # Switching Function Materials that provide h(eta_i)
    etas     = 'gr0 gr1' # Order parameters for all phases.
    # args = 'c' # Arguments of the fi free energies - use vector coupling
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
  [./top_displacement]
    type = FunctionDirichletBC
    variable = disp_y
    boundary = top
    function = ${my_function}
  [../]
  [./x_anchor]
    type = DirichletBC
    variable = disp_x
    boundary = 'left'
    value = 0.0
  [../]
  [./y_anchor]
    type = DirichletBC
    variable = disp_y
    boundary = 'bottom'
    value = 0.0
  [../]
[]

[Preconditioning]
  [./SMP]
    type = SMP
    full = true
  [../]
[]

# We monitor the total free energy and the total solute concentration (should be constant)
[Postprocessors]
#   [./total_free_energy]
#     type = ElementIntegralVariablePostprocessor
#     variable = local_energy
#   [../]
  [./dt]
    type = TimestepSize
  [../]
  [./active_time]     # Time computer spent on simulation
    type = PerfGraphData
    section_name = "Root"
    data_type = total
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
  end_time = 100

  [./TimeStepper]
    type = SolutionTimeAdaptiveDT
    dt = 0.1
  [../]
[]

[Outputs]
  file_base = ./${my_filename}/out_${my_filename} 
  execute_on = 'timestep_end'
  [./my_exodus]
    type = Exodus
    interval = 2
    # append_date = true
    # append_date_format = '%Y-%m-%d'
    # discontinuous = true
    sequence = true
  [../] 
  csv = true
  [./my_console]
    type = Console
    output_linear = false
    # output_screen = false
    # interval = 5
  [../]
  [./pgraph]
    type = PerfGraphOutput
    execute_on = 'initial timestep_end'  # Default is "final"
    level = 2                     # Default is 1
    heaviest_branch = true        # Default is false
    heaviest_sections = 7         # Default is 0
    # interval = 5
  [../]
[]

[Debug]
  # show_var_residual_norms = true
[]
