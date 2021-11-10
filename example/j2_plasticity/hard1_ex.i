# UserObject J2 test, with hardening, but with rate=0
# apply uniform compression in x direction to give
# trial stress_xx = -5, so sqrt(3*J2) = 5
# with zero Poisson's ratio, this should return to
# stress_xx = -3, stress_yy = -1 = stress_zz,
# for strength = 2
# (note that stress_xx - stress_yy = stress_xx - stress_zz = -2, so sqrt(3*j2) = 2,
#  and that the mean stress remains = -5/3)

my_num_element = 1
my_file_name = 'PowerRule'
[Mesh]
  displacements = 'disp_x disp_y'
  [./generated_mesh]
    type = GeneratedMeshGenerator
    dim = 2
    nx = ${my_num_element}
    ny = ${my_num_element}
    nz = ${my_num_element}
    xmin = 0
    xmax = 10
    ymin = 0
    ymax = 10
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

[Kernels]
  [./TensorMechanics]
    displacements = 'disp_x disp_y'
  [../]
[]


[BCs]
  [./x]
    type = DirichletBC
    variable = disp_x
    boundary = left
    value = 0
  [../]
  [./y]
    type = DirichletBC
    variable = disp_y
    boundary = bottom
    value = 0
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
    function = '0.01*t'
  [../]
[]

[AuxVariables]
#   [./strain_xx]
#     order = CONSTANT
#     family = MONOMIAL
#   [../]
#   [./strain_xy]
#     order = CONSTANT
#     family = MONOMIAL
#   [../]
  [./strain_yy]
    order = CONSTANT
    family = MONOMIAL
  [../]
#   [./plastic_strain_yy]
#     order = CONSTANT
#     family = MONOMIAL
#   [../]
#   [./stress_xx]
#     order = CONSTANT
#     family = MONOMIAL
#   [../]
#   [./stress_xy]
#     order = CONSTANT
#     family = MONOMIAL
#   [../]
  [./stress_yy]
    order = CONSTANT
    family = MONOMIAL
  [../]
#   [./VMstrain]
#     order = CONSTANT
#     family = MONOMIAL
#   [../]
  [./VMstress]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./f]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./p_internal_parameter]
    order = CONSTANT
    family = MONOMIAL
  [../]
#   [./iter]
#     order = CONSTANT
#     family = MONOMIAL
#   [../]
[]

[AuxKernels]
#   # active = 'none'
#   [./strain_xx]
#     type = RankTwoAux
#     rank_two_tensor = total_strain
#     variable = strain_xx
#     index_i = 0
#     index_j = 0
#   [../]
#   [./strain_xy]
#     type = RankTwoAux
#     rank_two_tensor = total_strain
#     variable = strain_xy
#     index_i = 0
#     index_j = 1
#   [../]
  [./strain_yy]
    type = RankTwoAux
    rank_two_tensor = total_strain
    variable = strain_yy
    index_i = 1
    index_j = 1
  [../]
#   [./plastic_strain_yy]
#     type = RankTwoAux
#     rank_two_tensor = plastic_strain
#     variable = plastic_strain_yy
#     index_i = 1
#     index_j = 1
#   [../]
#   [./stress_xx]
#     type = RankTwoAux
#     rank_two_tensor = stress
#     variable = stress_xx
#     index_i = 0
#     index_j = 0
#   [../]
#   [./stress_xy]
#     type = RankTwoAux
#     rank_two_tensor = stress
#     variable = stress_xy
#     index_i = 0
#     index_j = 1
#   [../]
  [./stress_yy]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_yy
    index_i = 1
    index_j = 1
  [../]
#   [./VMstrain]
#     type = MaterialRealAux
#     variable = VMstrain
#     property = equivalent_plastic_strain
#   [../]
  [./VMstress]
    type = MaterialRealAux
    variable = VMstress
    property = von_mises_stress  
  [../]
  [./f]
    type = MaterialStdVectorAux
    index = 0
    property = plastic_yield_function
    variable = f
  [../]
  [./p_internal_parameter]
    type = MaterialStdVectorAux
    index = 0
    property = plastic_internal_parameter
    variable = p_internal_parameter
  [../]
#   [./iter]
#     type = MaterialRealAux
#     property = plastic_NR_iterations
#     variable = iter
#   [../]
[]

[Postprocessors]
#   [./epsilon_xx]
#     type = PointValue
#     point = '0 0 0'
#     variable = strain_xx
#   [../]
#   [./epsilon_xy]
#     type = PointValue
#     point = '0 0 0'
#     variable = strain_xy
#   [../]
  [./epsilon_yy]
    type = PointValue
    point = '0 0 0'
    variable = strain_yy
  [../]
#   [./epsilon_p_yy]
#     type = PointValue
#     point = '0 0 0'
#     variable = plastic_strain_yy
#   [../]
#   [./sigma_xx]
#     type = PointValue
#     point = '0 0 0'
#     variable = stress_xx
#   [../]
#   [./sigma_xy]
#     type = PointValue
#     point = '0 0 0'
#     variable = stress_xy
#   [../]
  [./sigma_yy]
    type = PointValue
    point = '0 0 0'
    variable = stress_yy
  [../]
#   [./VMstrain]
#     type = PointValue
#     point = '0 0 0'
#     variable = VMstrain
#   [../]
  [./VMstress]
    type = PointValue
    point = '0 0 0'
    variable = VMstress
  [../]
  [./f]
    type = PointValue
    point = '0 0 0'
    variable = f
  [../]
  [./p_internal_parameter]
    type = PointValue
    point = '0 0 0'
    variable = p_internal_parameter
  [../]
#   [./iter]
#     type = PointValue
#     point = '0 0 0'
#     variable = iter
#   [../]
[]

[UserObjects]
  # [./str]
  #   type = GGTensorMechanicsHardeningConstant
  #   value = 700
  # [../]
  [./str]
    type = TensorMechanicsHardeningLinear
    value_0 = 33.3 # MPa
    HardFactor = 2.0e3
    # exponent = 1e1
    # value_0 * (p / epsilon0 + 1)^{exponent})
  [../]
  [./j2]
    type = GGTensorMechanicsPlasticJ2
    yield_strength = str
    yield_function_tolerance = 1E-5
    internal_constraint_tolerance = 1E-9
    max_iterations = 10
  [../]
[]

[Materials]
  # [./elasticity_tensor]
  #   type = ComputeElasticityTensor
  #   block = 0
  #   fill_method = symmetric9
  #   C_ijkl = '1.27e5 0.708e5 0.708e5 1.27e5 0.708e5 1.27e5 0.7355e5 0.7355e5 0.7355e5'
  # [../]
  [./elasticity_tensor]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = 1.1e5
    poissons_ratio = 0.31
  [../]
  [./strain]
    type = ComputeFiniteStrain
    block = 0
    displacements = 'disp_x disp_y'
  [../]
  [./mc]
    type = GGComputeMultiPlasticityStress
    block = 0
    ep_plastic_tolerance = 1E-9
    plastic_models = j2
    debug_fspb = crash
    # outputs = my_exodus
    # output_properties = 'plastic_internal_parameter'
    # elastic_strain_ij, plastic_strain, stress_ij,
  [../]
[]


[Executioner]
  end_time = 5
  dt = 0.02
  type = Transient
[]


[Outputs]
  file_base = ./hard1_ex/hard6_ex
  # exodus = false
  [./my_exodus]
    type = Exodus
  [../]
  [./csv]
    type = CSV
  [../]
[]
