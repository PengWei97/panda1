# UserObject J2 test, with hardening, but with rate=0, 
# apply uniform compression in x direction to give
# trial stress_xx = -5, so sqrt(3*J2) = 5
# with zero Poisson's ratio, this should return to
# stress_xx = -3, stress_yy = -1 = stress_zz,
# for strength = 2
# (note that stress_xx - stress_yy = stress_xx - stress_zz = -2, so sqrt(3*j2) = 2,
#  and that the mean stress remains = -5/3)
# is prefectly hardening

[Mesh]
  displacements = 'x_disp y_disp z_disp'
  [generated_mesh]
    type = GeneratedMeshGenerator
    elem_type = HEX8
    dim = 3
    nx = 1
    ny = 1
    nz = 1
    xmin = 0.0
    xmax = 1.0
    ymin = 0.0
    ymax = 1.0
    zmin = 0.0
    zmax = 1.0
  []
  [cnode]
    type = ExtraNodesetGenerator
    coord = '0.0 0.0 0.0'
    new_boundary = 6
    input = generated_mesh
  []
  [snode]
    type = ExtraNodesetGenerator
    coord = '1.0 0.0 0.0'
    new_boundary = 7
    input = cnode
  []
[]


[Variables]
  [./x_disp]
    order = FIRST
    family = LAGRANGE
  [../]
  [./y_disp]
    order = FIRST
    family = LAGRANGE
  [../]
  [./z_disp]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[Kernels]
  [./TensorMechanics]
    displacements = 'x_disp y_disp z_disp'
    use_displaced_mesh = true
  [../]
[]

[Functions]
  [./topfunc]
    type = ParsedFunction
    value = 't'
  [../]
[]

[BCs]
  [./bottom3]
    type = DirichletBC
    variable = z_disp
    boundary = 0
    value = 0.0
  [../]
  [./top]
    type = FunctionDirichletBC
    variable = z_disp
    boundary = 5
    function = topfunc
  [../]
  [./corner1]
    type = DirichletBC
    variable = x_disp
    boundary = 6
    value = 0.0
  [../]
  [./corner2]
    type = DirichletBC
    variable = y_disp
    boundary = 6
    value = 0.0
  [../]
  [./corner3]
    type = DirichletBC
    variable = z_disp
    boundary = 6
    value = 0.0
  [../]
  [./side1]
    type = DirichletBC
    variable = y_disp
    boundary = 7
    value = 0.0
  [../]
  [./side2]
    type = DirichletBC
    variable = z_disp
    boundary = 7
    value = 0.0
  [../]
[]

[AuxVariables]
#  active = ' '
  [./stress_zz]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./f]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./strain_zz]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./plastic_strain_zz]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./iter]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[AuxKernels]
#  active = ' '
  [./stress_zz]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_zz
    index_i = 2
    index_j = 2
  [../]
  [./strain_zz]
    type = RankTwoAux
    rank_two_tensor = total_strain
    variable = strain_zz
    index_i = 2
    index_j = 2
  [../]
  [./plastic_strain_zz]
    type = RankTwoAux
    rank_two_tensor = plastic_strain
    variable = plastic_strain_zz
    index_i = 2
    index_j = 2
  [../]
  [./f]
    type = MaterialStdVectorAux
    index = 0
    property = plastic_yield_function
    variable = f
  [../]
  [./iter]
    type = MaterialRealAux
    property = plastic_NR_iterations
    variable = iter
  [../]
[]

[Postprocessors]
  # active = ' '
  [./total_stress_zz]
    type = ElementAverageValue
    variable = stress_zz
  [../]
  [./total_strain_zz]
    type = ElementAverageValue
    variable = strain_zz
  [../]
  [./f]
    type = ElementAverageValue
    variable = f
  [../]
  [./iter]
    type = PointValue
    point = '0 0 0'
    variable = iter
  [../]
[]

[UserObjects]
  [./str]
    type = TensorMechanicsHardeningConstant
    value = 4e5
  [../]
  [./j2]
    type = TensorMechanicsPlasticJ2
    yield_strength = str
    yield_function_tolerance = 1E-5
    internal_constraint_tolerance = 1E-9
  [../]
[]

[Materials]
  [./elasticity_tensor]
    type = ComputeElasticityTensor
    block = 0
    fill_method = symmetric_isotropic
    C_ijkl = '0 1E6'
    # elasticity_tensor_ijkl, effective_stiffness
  [../]
  [./strain]
    type = ComputeFiniteStrain
    block = 0
    displacements = 'x_disp y_disp z_disp' 
    # deformation_gradient_ij, mechanical_strain_ij
    # rotation_increment_ij, strain_increment_ij
    # total_strain_ij.strain_rate_ij
  [../]
  [./mc]
    type = ComputeMultiPlasticityStress
    block = 0
    ep_plastic_tolerance = 1E-9
    plastic_models = j2
    debug_fspb = crash
    # outputs = exodus
    # Jacobian_mult_ijkl, elastic_strain_ij,
    # plastic_NR_iterations
    # plastic_strain_ij
    # stress_ij
  [../]
[]

[Preconditioning]
  [./SMP]
   type = SMP
   full=true
  [../]
[]

[Executioner]
  type = Transient

  dt=0.01
  # dtmax=1
  dtmin=0.01
  end_time=2

  nl_abs_tol = 1e-10
[]


[Outputs]
  file_base = ./j2_constHardening/j2_const6_out
  exodus = true
  [./csv]
    type = CSV
  [../]
[]
