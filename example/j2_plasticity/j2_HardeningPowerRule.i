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
  displacements = 'x_disp y_disp' # z_disp
  [generated_mesh]
    type = GeneratedMeshGenerator
    elem_type = QUAD4
    dim = 2
    nx = 20
    ny = 20
    nz = 0
    xmin = 0.0
    xmax = 10.0
    ymin = 0.0
    ymax = 10.0
    zmin = 0.0
    zmax = 0.0
  []
  [cnode]
    type = ExtraNodesetGenerator
    coord = '0.0 0.0'
    new_boundary = 6
    input = generated_mesh
  []
  [snode]
    type = ExtraNodesetGenerator
    coord = '1.0 0.0'
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
  # [./z_disp]
  #   order = FIRST
  #   family = LAGRANGE
  # [../]
[]

[Kernels]
  [./TensorMechanics]
    displacements = 'x_disp y_disp'
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
    variable = y_disp
    boundary = bottom
    value = 0.0
  [../]
  [./top]
    type = FunctionDirichletBC
    variable = y_disp
    boundary = top
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
  # [./corner3]
  #   type = DirichletBC
  #   variable = z_disp
  #   boundary = 6
  #   value = 0.0
  # [../]
  [./right2]
    type = DirichletBC
    variable = x_disp
    boundary = right
    value = 0.0
  [../]
  # [./side2]
  #   type = DirichletBC
  #   variable = z_disp
  #   boundary = 7
  #   value = 0.0
  # [../]
[]

[AuxVariables]
#  active = ' '
  [./stress_yy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./f]
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
  [./iter]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[AuxKernels]
#  active = ' '
  [./stress_yy]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_yy
    index_i = 1
    index_j = 1
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
  [./total_stress_yy]
    type = ElementAverageValue
    variable = stress_yy
  [../]
  [./total_strain_yy]
    type = ElementAverageValue
    variable = strain_yy
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
    type = TensorMechanicsHardeningPowerRule
    value_0 = 4e5
    epsilon0 = 0.1
    exponent = 0.1
  [../]
  [./j2]
    type = TensorMechanicsPlasticJ2
    yield_strength = str
    yield_function_tolerance = 1E-5
    internal_constraint_tolerance = 1E-9
  [../]
[]

# value = value_0 * (p / epsilon0 + 1)^{exponent})
# sigma_y_0 = 4e5;
# epsilon0 = 0.1;
# n = 0.1
# eq_epsilon = linspace(0,2,100)
# sigma_y = sigma_y_0.*(eq_epsilon./epsilon0+1).^n
# plot(eq_epsilon,sigma_y)

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
    displacements = 'x_disp y_disp' 
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
  end_time=3

  nl_abs_tol = 1e-10
  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.01
    growth_factor = 1.2
    cutback_factor = 0.8
    optimal_iterations = 8
  [../]
  # [./Adaptivity]
  #   initial_adaptivity = 3 # 8 
  #   cycles_per_step = 2 # The number of adaptivity cycles per step
  #   refine_fraction = 0.5 # The fraction of elements or error to refine.
  #   coarsen_fraction = 0.05
  #   max_h_level = 8
  # [../]
[]



[Outputs]
  file_base = ./j2_HardeningPowerRule/j2_HardeningPowerRule_2d_out
  exodus = true
  [./csv]
    type = CSV
  [../]
[]
