my_filename = 'test03_phJ2_01'
# my_function = '100+0.01*sin(t*pi)'
# my_function = '10*t' # 0.001s^{-1}
my_function = 'if(t<40,0.1*t,40+0.02*sin(t))' # 0.001s^{-1}
# my_function = 'if(t<9,20*t,180+0.1*sin(t*pi))' # 0.001s^{-1}
# my_function = 'if(t<3,20*t,60+0.02*sin(t*pi))'
my_end_time = 1e5
# my_yield_0 = 700 # MPa
my_xmax = 3e3
my_ymax = 1e3
my_radus = 1e3

my_time_scale = 1.0e-9
my_length_scale = 1.0e-9
my_pressure_scale = 1.0e6

my_yield_strength_init = 2e3

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 50
  ny = 20
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

[Variables]
  [./PolycrystalVariables]
  [../]
  [./disp_x]
  [../]
  [./disp_y]
  [../]
[]

[ICs]
  [./PolycrystalICs]
    [./BicrystalBoundingBoxIC]
      x1 = 0
      y1 = 0
      x2 = ${my_radus}
      y2 = ${my_ymax}   
    [../]
  [../]
[]

[AuxVariables]
  [./bnds]
    order = FIRST
    family = LAGRANGE
  [../]
  [./strain_xx]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./strain_yy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./elastic_strain11]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./elastic_strain22]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_xx]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_yy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./VMstress]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./p_internal_parameter]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./HardFactor]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./unique_grains]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./var_indices]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./C1111]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./euler_angle]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[Kernels]
  [./PolycrystalKernel]
  [../]
  [./GGPolycrystalElasticDrivingForce]
    # ACGrGrElasticDrivingForce
    # GGPolycrystalElasticDrivingForce
    op_num = 2
    var_name_base = gr
      # GGACGrGrElasticDrivingForce
  [../]
  # [./ACGrGrPlasticDrivingForce1]
  #   type = ACGrGrPlasticDrivingForce
  #   D_hard_factor_name = dhard_factor/dgr0
  #   variable = gr0
  # [../]
  # [./ACGrGrPlasticDrivingForce2]
  #   type = ACGrGrPlasticDrivingForce
  #   D_hard_factor_name = dhard_factor/dgr1
  #   variable = gr1
  # [../]
[]

[Modules/TensorMechanics/Master]
  [./all]
    # strain = FINITE
    # displacements = 'disp_x disp_y'  
    # use_displaced_mesh = true
    strain = FINITE # FINITE
  [../]
[]

[AuxKernels]
  [./bnds_aux]
    type = BndsCalcAux
    variable = bnds
    execute_on = timestep_end
  [../]
  [./elastic_strain11]
    type = RankTwoAux
    variable = elastic_strain11
    rank_two_tensor = elastic_strain
    index_i = 0
    index_j = 0
    execute_on = timestep_end
  [../]
  [./elastic_strain22]
    type = RankTwoAux
    variable = elastic_strain22
    rank_two_tensor = elastic_strain
    index_i = 1
    index_j = 1
    execute_on = timestep_end
  [../]
  [./strain_xx]
    type = RankTwoAux
    rank_two_tensor = total_strain
    variable = strain_xx
    index_i = 0
    index_j = 0
  [../]
  [./strain_yy]
    type = RankTwoAux
    rank_two_tensor = total_strain
    variable = strain_yy
    index_i = 1
    index_j = 1
  [../]
  [./stress_xx]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_xx
    index_i = 0
    index_j = 0
  [../]
  [./stress_yy]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_yy
    index_i = 1
    index_j = 1
  [../]
  [./VMstress]
    type = MaterialRealAux
    variable = VMstress
    property = von_mises_stress  
  [../]
  [./p_internal_parameter]
    type = MaterialRealAux
    # index = 0
    property = eqv_plastic_strain
    variable = p_internal_parameter
  [../]
  [./HardFactor]
    type = MaterialRealAux
    variable = HardFactor
    property = hard_factor  
  [../]
  [./unique_grains]
    type = FeatureFloodCountAux
    variable = unique_grains
    flood_counter = grain_tracker
    execute_on = 'initial timestep_begin'
    field_display = UNIQUE_REGION
  [../]
  [./var_indices]
    type = FeatureFloodCountAux
    variable = var_indices
    flood_counter = grain_tracker
    execute_on = 'initial timestep_begin'
    field_display = VARIABLE_COLORING
  [../]
  [./C1111]
    type = RankFourAux
    variable = C1111
    rank_four_tensor = elasticity_tensor
    index_l = 0
    index_j = 0
    index_k = 0
    index_i = 0
    execute_on = timestep_end
  [../]
  [./euler_angle]
    type = OutputEulerAngles
    variable = euler_angle
    euler_angle_provider = euler_angle_file
    grain_tracker = grain_tracker
    output_euler_angle = 'phi1'
  [../]
[]

[BCs]
  [./top_displacement]
    type = FunctionDirichletBC
    variable = disp_y
    boundary = top
    function = ${my_function}
    use_displaced_mesh = true
    # function = 50
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

[Materials]
  [./Copper]
    type = GBEvolution
    block = 0
    T = 500 # K
    wGB = 75 # nm
    GBmob0 = 2.5e-6 #m^4/(Js) from Schoenfelder 1997
    Q = 0.23 # Migration energy in eV
    GBenergy = 0.708 # GB energy in J/m^2
    time_scale = ${my_time_scale}
    length_scale = ${my_length_scale}
  [../]
  [./ElasticityTensor]
    type = GGComputePolycrystalElasticityTensor
    grain_tracker = grain_tracker
  [../]
  [./fplastic]
    type = Test4FiniteStrainPlasticMaterial #设置屈服函数
    # implements rate-independent associative J2 plasticity 
    # with isotropic hardening in the finite-strain framework.
    block = 0
    grain_tracker = grain_tracker
    yield_stress='0. 700. 0.05 700. 0.1 700. 0.38 700. 0.95 700. 2. 700.'
    outputs = my_exodus
    output_properties = 'dhard_factor/dgr0 dhard_factor/dgr1'
    length_scale = ${my_length_scale}
    pressure_scale = ${my_pressure_scale}
    yield_strength_init = ${my_yield_strength_init}
  [../]
  [./elastic_free_energy]
    type = ElasticEnergyMaterial
    f_name = f_elastic
    block = 0
    args = 'gr0 gr1'
    outputs = my_exodus
    output_properties = 'f_elastic df_elastic/dgr0 df_elastic/dgr1'
  [../]
  [./local_free_energy]
    type = DerivativeParsedMaterial
    f_name= f_chem
    args = 'gr0 gr1'
    material_property_names = 'mu gamma_asymm'
    function = 'mu*(gr0^4/4.0 - gr0^2/2.0 + gr1^4/4.0 - gr1^2/2.0 + gamma_asymm*gr0^2*gr1^2+1.0/4.0)'
    derivative_order = 2
    enable_jit = true
    outputs = my_exodus
    output_properties = 'f_chem df_chem/dgr0 df_chem/dgr1'
  [../]
[]

[UserObjects]
  [./euler_angle_file]
    type = EulerAngleFileReader
    file_name = test.tex
  [../]
  [./grain_tracker]
    type = GrainTrackerElasticity
    connecting_threshold = 0.05
    compute_var_to_feature_map = true
    flood_entity_type = elemental
    execute_on = 'initial timestep_begin'

    euler_angle_provider = euler_angle_file
    fill_method = symmetric9
    C_ijkl = '1.684e5 1.214e5 1.214e5 1.684e5 1.214e5 1.684e5 0.75e5 0.75e5 0.75e5'

    outputs = none
  [../]
[]

[Postprocessors]
  [./dt]
    type = TimestepSize
  [../]
  [./gr0_area]
    type = ElementIntegralVariablePostprocessor
    variable = gr0
  [../]
  [./epsilon_yy]
    type = PointValue
    point = '0 0 0'
    variable = strain_yy
  [../]
  [./sigma_yy]
    type = PointValue
    point = '0 0 0'
    variable = stress_yy
  [../]
  [./VMstress]
    type = PointValue
    point = '0 0 0'
    variable = VMstress
  [../]
  [./p_internal_parameter]
    type = PointValue
    point = '0 0 0'
    variable = p_internal_parameter
  [../]
  [./HardFactor]
    type = PointValue
    point = '0 0 0'
    variable = HardFactor
  [../]
  # [./active_time]           # Time computer spent on simulation
  #   type = PerfGraphData
  #   section_name = "Root"
  #   data_type = total
  # [../]
[]

# [Preconditioning]
#   [./SMP]
#    type = SMP
#    coupled_groups = 'gr0,gr1 disp_x,disp_y'
#   [../]
# []

[Preconditioning]
  [./SMP]
   type = SMP
   coupled_groups = 'disp_x,disp_y'
  [../]
[]

[Executioner]
  type = Transient

  solve_type = 'PJFNK'
  # petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart -pc_hypre_boomeramg_strong_threshold'
  # petsc_options_value = 'hypre boomeramg 31 0.7'

  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'

  l_max_its = 30
  l_tol = 1e-4
  nl_max_its = 30
  nl_rel_tol = 1e-9

  start_time = 0.0
  end_time = ${my_end_time}
  # # num_steps = 30
  dt = 0.1
  # [./TimeStepper]
  #   type = IterationAdaptiveDT
  #   dt = 0.2
  #   growth_factor = 1.2
  #   cutback_factor = 0.8
  #   optimal_iterations = 8
  # [../]
  # [./Adaptivity]
  #   initial_adaptivity = 5
  #   refine_fraction = 0.7
  #   coarsen_fraction = 0.1
  #   max_h_level = 3
  # [../]
[]

[Outputs]
  file_base = ./${my_filename}/out_${my_filename} 
  execute_on = 'timestep_end'
  [./my_exodus]
    type = Exodus
  [../] 
  csv = true
  [./my_console]
    type = Console
    output_linear = false
    # output_screen = false
    interval = 5
  [../]
  [./pgraph]
    type = PerfGraphOutput
    execute_on = 'initial timestep_end final'  # Default is "final"
    level = 2                     # Default is 1
    heaviest_branch = true        # Default is false
    heaviest_sections = 7         # Default is 0
  [../]
[]
