my_filename = 'test03_j2Ph_circular_5'

# test03_j2Ph_circular_3 外面晶粒为45°，环晶长大
# test03_j2Ph_circular_4 里面晶粒为45°，环晶长大, 晶粒缩小不明显
# test03_j2Ph_circular_5 里面晶粒为45°，环晶长大，设置载荷为0
my_xnum_element = 32
my_ynum_element = 32
my_xmax = 6.4e2
my_ymax = 6.4e2
# my_function = 't' # 0.001
# my_function = 'if(t<40,0.64*t,25.6+0.02*sin(t))' # 0.001s^{-1} 4%
my_function = '0' # 0.001s^{-1} 4%
my_end_time = 1.0e3
my_HardFactor = 2e4 # 30000.0

my_radus = 2.5e2
my_point = 3.2e2
my_time_scale = 1.0e-9
my_length_scale = 1.0e-9
my_pressure_scale = 1.0e6

my_num_adaptivity = 3
[Mesh]
  displacements = 'disp_x disp_y'
  [./generated_mesh]
    type = GeneratedMeshGenerator
    dim = 2
    nx = ${my_xnum_element}
    ny = ${my_ynum_element}
    xmax = ${my_xmax}
    ymax = ${my_ymax}
    elem_type = QUAD4
    block = 0
  [../]
  [cnode]
    type = ExtraNodesetGenerator
    coord = '0.0 0.0'
    new_boundary = 6
    input = generated_mesh
  []
[]

[GlobalParams]
  op_num = 2
  var_name_base = gr
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
    [./BicrystalCircleGrainIC]
      radius = ${my_radus}
      x = ${my_point}
      y = ${my_point}
      int_width = 15
    [../]
  [../]
[]

[BCs]
  [./x]
    type = DirichletBC
    variable = disp_x
    boundary = 'left'
    value = 0.0
  [../]
  [./y]
    type = DirichletBC
    variable = disp_y
    boundary = 'bottom'
    value = 0.0
  [../]
  [./corner1]
    type = DirichletBC
    variable = disp_x
    boundary = 6
    value = 0.0
  [../]
  [./corner2]
    type = DirichletBC
    variable = disp_y
    boundary = 6
    value = 0.0
  [../]
  [./top]
    type = FunctionDirichletBC
    variable = disp_y
    boundary = top
    function = ${my_function}
  [../]
[]

[AuxVariables]
  [./bnds]
    order = FIRST
    family = LAGRANGE
  [../]
  [./strain_yy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./plastic_strain_yy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./elastic_strain_yy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./elastic_strain_xx]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_yy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_xx]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./VMstress]
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

[AuxKernels]
  [./strain_yy]
    type = RankTwoAux
    rank_two_tensor = total_strain
    variable = strain_yy
    index_i = 1
    index_j = 1
  [../]
  [./plastic_strain_yy]
    type = RankTwoAux
    rank_two_tensor = plastic_strain
    variable = plastic_strain_yy
    index_i = 1
    index_j = 1
  [../]
  [./elastic_strain_yy]
    type = RankTwoAux
    rank_two_tensor = elastic_strain
    variable = elastic_strain_yy
    index_i = 1
    index_j = 1
  [../]
  [./elastic_strain_xx]
    type = RankTwoAux
    rank_two_tensor = elastic_strain
    variable = elastic_strain_xx
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
  [./stress_xx]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_xx
    index_i = 0
    index_j = 0
  [../]
  [./VMstress]
    type = MaterialRealAux
    variable = VMstress
    property = von_mises_stress  
  [../]
  [./bnds_aux]
    type = BndsCalcAux
    variable = bnds
    execute_on = timestep_end
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



[Postprocessors]
  [./epsilon_yy]
    type = PointValue
    point = '0 0 0'
    variable = strain_yy
  [../]
  [./epsilon_yy_av]
    type = ElementAverageValue
    variable = strain_yy
  [../]
  [./epsilon_p_yy]
    type = PointValue
    point = '0 0 0'
    variable = plastic_strain_yy
  [../]
  [./epsilon_e_yy]
    type = PointValue
    point = '0 0 0'
    variable = elastic_strain_yy
  [../]
  [./epsilon_e_xx]
    type = PointValue
    point = '0 0 0'
    variable = elastic_strain_xx
  [../]
  [./sigma_yy]
    type = PointValue
    point = '0 0 0'
    variable = stress_yy
  [../]
  [./sigma_xx]
    type = PointValue
    point = '0 0 0'
    variable = stress_xx
  [../]
  [./VMstress]
    type = PointValue
    point = '0 0 0'
    variable = VMstress
  [../]
  [./VMstress_av]
    type = ElementAverageValue
    variable = VMstress
  [../]
  [./run_time]
    type = PerfGraphData
    section_name = "Root"
    data_type = total
  [../]
  [./dt]
    type = TimestepSize
  [../]
  [./gr0_area]
    type = ElementIntegralVariablePostprocessor
    variable = gr0
  [../]
  [./F_elastic]
    type = ElementIntegralMaterialProperty
    mat_prop = f_elastic
  [../]
  [./F_chem]
    type = ElementIntegralMaterialProperty
    mat_prop = f_chem
    # outputs = csv
  [../]
[]

[UserObjects]
  [./str]
    type = TensorMechanicsHardeningLinear
    value_0 = 2000 # MPa
    HardFactor = ${my_HardFactor} # 30000.0 # MPa
  [../]
  [./j2]
    type = GGTensorMechanicsPlasticJ2
    yield_strength = str
    yield_function_tolerance = 1E-5
    internal_constraint_tolerance = 1E-9
    max_iterations = 10
  [../]
  [./euler_angle_file]
    type = EulerAngleFileReader
    file_name = Euler_circular.tex
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

[Modules/TensorMechanics/Master]
  [./all]
    # strain = FINITE     
    # use_displaced_mesh = true
    strain = FINITE # FINITE
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
[]

[Materials]
  [./mc]
    type = GGComputeMultiPlasticityStress
    block = 0
    ep_plastic_tolerance = 1E-9
    plastic_models = j2
    debug_fspb = crash
    # tangent_operator = elastic
    # perform_finite_strain_rotations = false
    # plastic_strain, plastic_internal_parameter, 
    # plastic_yield_function
    # total_strain, stress, elastic_strain
  [../]
  [./Copper]
    type = GBEvolution
    block = 0
    T = 500 # K
    wGB = 75 # nm
    GBmob0 = 2.5e-6 #m^4/(Js) from Schoenfelder 1997
    Q = 0.23 #Migration energy in eV
    GBenergy = 0.708 #GB energy in J/m^2
    time_scale = ${my_time_scale}
    length_scale = ${my_length_scale}
  [../]
  # [./HardeningModel]
  #   type = GGTensorMechanicsHardeningLinear
  #   # grain_tracker = grain_tracker
  # [../]
  [./ElasticityTensor]
    type = GGComputePolycrystalElasticityTensor
    grain_tracker = grain_tracker
  [../]
  [./elastic_free_energy]
    type = GGElasticEnergyMaterial
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

[Preconditioning]
  [./SMP]
   type = SMP
   coupled_groups = 'gr0 gr1 disp_x,disp_y'
  [../]
[]

[Executioner]
  type = Transient

  solve_type = 'PJFNK'
  petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart -pc_hypre_boomeramg_strong_threshold'
  petsc_options_value = 'hypre boomeramg 31 0.7'
  # petsc_options_iname = '-pc_type' # Names of PETSc name/value pairs
  # petsc_options_value = 'lu' # 

  dtmin = 2.0e-3 # The minimum timestep size in an adaptive run
  l_max_its = 30
  l_tol = 1e-4
  nl_max_its = 30
  nl_rel_tol = 1e-9

  end_time = ${my_end_time}
  dt = 0.2

  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.2
    growth_factor = 1.2
    cutback_factor = 0.8
    optimal_iterations = 8
  [../]

  [./Adaptivity]
    initial_adaptivity = ${my_num_adaptivity}
    refine_fraction = 0.7
    coarsen_fraction = 0.1
    max_h_level = ${my_num_adaptivity}
  [../]
[]


[Outputs]
  file_base = ./${my_filename}/out_${my_filename} 
  # file_base = ./test005/out_${my_filename}
  [./my_exodus]
    type = Exodus
  [../]
  [./csv]
    type = CSV
  [../]
  [./pgraph]
    type = PerfGraphOutput
    execute_on = 'initial timestep_end final'  # Default is "final"
    level = 2                     # Default is 1
    heaviest_branch = true        # Default is false
    heaviest_sections = 7         # Default is 0
  [../]
[]
