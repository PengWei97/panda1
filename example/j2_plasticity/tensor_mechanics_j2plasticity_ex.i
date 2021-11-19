# UserObject J2 test, with hardening, but with rate=0
# apply uniform compression in x direction to give
# trial stress_xx = -5, so sqrt(3*J2) = 5
# with zero Poisson's ratio, this should return to
# stress_xx = -3, stress_yy = -1 = stress_zz,
# for strength = 2
# (note that stress_xx - stress_yy = stress_xx - stress_zz = -2, so sqrt(3*j2) = 2,
#  and that the mean stress remains = -5/3)
my_filename = 'test06_2'
my_xnum_element = 60
my_ynum_element = 60
my_xmax = 3e3
my_ymax = 1e3
# my_function = 'if(t<6,10*t,60+0.01*sin(t))'
# my_function = 'if(t<3,20*t,60+0.02*sin(t*pi))'
my_function = 't' # 0.001
my_end_time = 150

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
  [./disp_x]
  [../]
  [./disp_y]
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
  [./VMstrain]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./VMstress]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./HardFactor]
    order = CONSTANT
    family = MONOMIAL
  [../]
  # [./f]
  #   order = CONSTANT
  #   family = MONOMIAL
  # [../]
  # [./iter]
  #   order = CONSTANT
  #   family = MONOMIAL
  # [../]
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
  [./stress_yy]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_yy
    index_i = 1
    index_j = 1
  [../]
  [./VMstrain]
    type = MaterialRealAux
    variable = VMstrain
    property = eqv_plastic_strain
  [../]
  [./VMstress]
    type = MaterialRealAux
    variable = VMstress
    property = von_mises_stress  
  [../]
  [./HardFactor]
    type = MaterialRealAux
    variable = HardFactor
    property = hard_factor  
  [../]
  # [./f]
  #   type = MaterialStdVectorAux
  #   index = 0
  #   property = plastic_yield_function
  #   variable = f
  # [../]
  # [./iter]
  #   type = MaterialRealAux
  #   property = plastic_NR_iterations
  #   variable = iter
  # [../]
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
  [./sigma_yy]
    type = PointValue
    point = '0 0 0'
    variable = stress_yy
  [../]
  [./VMstrain]
    type = PointValue
    point = '0 0 0'
    variable = VMstrain
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
  [./HardFactor]
    type = PointValue
    point = '0 0 0'
    variable = HardFactor
  [../]
  [./run_time]
    type = PerfGraphData
    section_name = "Root"
    data_type = total
  [../]
  [./dt]
    type = TimestepSize
  [../]
  # [./f]
  #   type = PointValue
  #   point = '0 0 0'
  #   variable = f
  # [../]
  # [./iter]
  #   type = PointValue
  #   point = '0 0 0'
  #   variable = iter
  # [../]
[]

# [UserObjects]
#   [./str]
#     type = GGTensorMechanicsHardeningConstant
#     value = 700
#   [../]
#   # [./str]
#   #   type = TensorMechanicsHardeningPowerRule
#   #   value_0 = 700
#   #   epsilon0 = 1
#   #   exponent = 1e1
#   #   # value_0 * (p / epsilon0 + 1)^{exponent})
#   # [../]
#   [./j2]
#     type = GGTensorMechanicsPlasticJ2
#     yield_strength = str
#     yield_function_tolerance = 1E-5
#     internal_constraint_tolerance = 1E-9
#     max_iterations = 10
#   [../]
# []

[Modules/TensorMechanics/Master]
  [./all]
    # strain = FINITE     
    # use_displaced_mesh = true
    strain = FINITE # FINITE
  [../]
[]

[GlobalParams]
  displacements = 'disp_x disp_y'
[]

[Materials]
  [./elasticity_tensor]
    type = ComputeElasticityTensor
    block = 0
    fill_method = symmetric9
    C_ijkl = '1.27e5 0.708e5 0.708e5 1.27e5 0.708e5 1.27e5 0.7355e5 0.7355e5 0.7355e5'
    # use_displaced_mesh = true
  [../]
  [./fplastic]
    type = Test2FiniteStrainPlasticMaterial # 设置屈服函数
    # implements rate-independent associative J2 plasticity 
    # with isotropic hardening in the finite-strain framework.
    block = 0
    yield_stress='0. 445. 0.05 610. 0.1 680. 0.38 810. 0.95 920. 2. 950.'
    # use_displaced_mesh = true
    # outputs = my_exodus
    # output_properties = 'hard_factor'
    # (610-445)/0.05 = 3300
    # value_0 = 445
    # hard_factor = 3300
    # equivalent plastic strain = ??
    # Yield function = sqrt(3*s_ij*s_ij/2) - K(equivalent plastic strain)
    # s_ij = stress_ij - delta_ij*trace(stress)/3
    # declareProperty
      # plastic_strain
      # eqv_plastic_strain
    # getMaterialProperty(old)
      # plastic_strain
      # eqv_plastic_strain
      # stress
      # strain_increment
      # rotation_increment
      # elasticity_tensor    
  [../]
[]

[Executioner]
  type = Transient
  end_time = ${my_end_time}
  dt = 0.05
  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.2
    growth_factor = 1.2
    cutback_factor = 0.8
    optimal_iterations = 8
  [../]
[]


[Outputs]
  file_base = ./${my_filename}/out_${my_filename} 
  # exodus = false
  [./my_exodus]
    type = Exodus
  [../]
  [./csv]
    type = CSV
  [../]
[]
