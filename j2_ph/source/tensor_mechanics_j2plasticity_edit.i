my_filename = 'test02_j2_04'
my_xnum_element = 50
my_ynum_element = 20
my_xmax = 1 #3e3
my_ymax = 1 #1e3

# my_function = 't' # 0.001
my_function = 'if(t<40,t,40+0.02*sin(t))' # 0.001s^{-1}
my_end_time = 1000
my_HardFactor = 2e4 # 30000.0

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
    order = FIRST
    family = LAGRANGE
  [../]
  [./disp_y]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[Kernels]
  [./TensorMechanics]
    displacements = 'disp_x disp_y'
    # use_displaced_mesh = true
  [../]
[]

[Materials]
  [./fplastic]
    type = FiniteStrainPlasticMaterial #设置屈服函数
    # implements rate-independent associative J2 plasticity 
    # with isotropic hardening in the finite-strain framework.
    block=0
    yield_stress='0. 445. 0.05 610. 0.1 680. 0.38 810. 0.95 920. 2. 950.'
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
  [./elasticity_tensor]
    type = ComputeElasticityTensor
    block = 0
    C_ijkl = '2.827e5 1.21e5 1.21e5 2.827e5 1.21e5 2.827e5 0.808e5 0.808e5 0.808e5'
    fill_method = symmetric9
  [../]
  [./strain]
    type = ComputeFiniteStrain
    block = 0
    displacements = 'disp_x disp_y'
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


# [AuxVariables]
#   [./stress_zz]
#     order = CONSTANT
#     family = MONOMIAL
#   [../]
#   [./peeq]
#     order = CONSTANT
#     family = MONOMIAL
#   [../]
#   [./pe11]
#     order = CONSTANT
#     family = MONOMIAL
#   [../]
#   [./pe22]
#     order = CONSTANT
#     family = MONOMIAL
#   [../]
#   [./pe33]
#     order = CONSTANT
#     family = MONOMIAL
#   [../]
# []

# [AuxKernels]
#   [./stress_zz]
#     type = RankTwoAux
#     rank_two_tensor = stress
#     variable = stress_zz
#     index_i = 2
#     index_j = 2
#   [../]
#   [./pe11]
#     type = RankTwoAux
#     rank_two_tensor = plastic_strain
#     variable = pe11
#     index_i = 0
#     index_j = 0
#   [../]
#   [./pe22]
#     type = RankTwoAux
#     rank_two_tensor = plastic_strain
#     variable = pe22
#     index_i = 1
#     index_j = 1
#   [../]
#   [./pe33]
#     type = RankTwoAux
#     rank_two_tensor = plastic_strain
#     variable = pe33
#     index_i = 2
#     index_j = 2
#   [../]
#   [./eqv_plastic_strain]
#     type = MaterialRealAux
#     property = eqv_plastic_strain
#     variable = peeq
#   [../]
# []

# [Postprocessors]
#   [./eqv_plastic_strain]
#     type = PointValue
#     point = '0 0 0'
#     variable = peeq
#   [../]
# []

[Preconditioning]
  [./SMP]
   type = SMP
   full=true
  [../]
[]

[Executioner]
  type = Transient

  solve_type = 'PJFNK'
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'

  dtmin = 2.0e-6 # The minimum timestep size in an adaptive run
  l_max_its = 30
  l_tol = 1e-4
  nl_max_its = 30
  nl_rel_tol = 1e-9

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
  [./my_exodus]
    type = Exodus
  [../]
  [./csv]
    type = CSV
  [../]
[]
