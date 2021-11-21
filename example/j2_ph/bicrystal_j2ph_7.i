my_filename = 'test07_j2_03' # test05系列用于建立j2弹性能的耦合，主要设置硬化参数的梯度化
my_xnum_element = 50
my_ynum_element = 20
my_xmax = 3e3
my_ymax = 1e3
# my_function = 't' # 0.001
my_function = 'if(t<40,t,40+0.02*sin(t))' # 0.001s^{-1}
my_end_time = 1000
my_HardFactor = 2e4 # 30000.0

my_radus = 1.0e3

my_length_scale = 1.0e-9
my_time_scale = 1.0e-9

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
  [./euler_angle]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./C1111]
    order = CONSTANT
    family = MONOMIAL
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
  [./stress_yy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./VMstress]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[AuxKernels]
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
[]

[Functions]
  [./timestep_fn]
    type = ParsedFunction
    value = 'if(t<40,0.2,10)'
  [../]
[]

[Postprocessors]
  [./timestep_pp]
    type = FunctionValuePostprocessor
    function = timestep_fn
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
  [./VMstress_av]
    type = ElementAverageValue
    variable = VMstress
  [../]
  # [./run_time]
  #   type = PerfGraphData
  #   section_name = "Root"
  #   data_type = total
  # [../]
  [./dt]
    type = TimestepSize
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

[GlobalParams]
  op_num = 2
  var_name_base = gr
  displacements = 'disp_x disp_y'
[]

[Materials]
  # [./elasticity_tensor]
  #   type = ComputeElasticityTensor
  #   block = 0
  #   fill_method = symmetric9
  #   C_ijkl = '1.27e5 0.708e5 0.708e5 1.27e5 0.708e5 1.27e5 0.7355e5 0.7355e5 0.7355e5'
  #   # use_displaced_mesh = true
  # [../]
  [./fplastic]
    type = Test6FiniteStrainPlasticMaterial #设置屈服函数
    # implements rate-independent associative J2 plasticity 
    # with isotropic hardening in the finite-strain framework.
    block=0
    yield_stress='0. 445. 0.05 610. 0.1 680. 0.38 810. 0.95 920. 2. 950.'
    grain_tracker = grain_tracker
  [../]
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
[]

[Preconditioning]
  [./SMP]
   type = SMP
   coupled_groups = 'gr0,gr1 disp_x,disp_y'
  [../]
[]

[Executioner]
  type = Transient

  solve_type = 'PJFNK'
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'

  # petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart -pc_hypre_boomeramg_strong_threshold'
  # petsc_options_value = 'hypre boomeramg 31 0.7'

  dtmin = 2.0e-6 # The minimum timestep size in an adaptive run
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
    timestep_limiting_postprocessor = timestep_pp
    # time_t = '0 60 1000'
    # time_dt = '0.2 0.2 10'
  [../]
  [./Adaptivity]
    initial_adaptivity = 5
    refine_fraction = 0.7
    coarsen_fraction = 0.1
    max_h_level = 3
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
