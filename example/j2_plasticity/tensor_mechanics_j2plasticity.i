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

[Materials]
  [./fplastic]
    type = FiniteStrainPlasticMaterial # 设置屈服函数
    # implements rate-independent associative J2 plasticity 
    # with isotropic hardening in the finite-strain framework.
    block=0
    yield_stress='0. 445. 0.05 610. 0.1 680. 0.38 810. 0.95 920. 2. 950.'
    # equivalent plastic strain = \int_{0}^{t} \sqrt{\frac{2}{3} \dot{\varepsilon}^{p}: \dot{\varepsilon}^{p}} \mathrm{~d} t
    # Yield function = sqrt(3*s_ij*s_ij/2) - K(equivalent plastic strain)
    # s_ij = stress_ij - delta_ij*trace(stress)/3
    # K is the yield stress, specified as a piecewise-linear function by the user.
    # This is a piecewise linear function entered by the user in the yield_stress vector
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
    displacements = 'x_disp y_disp z_disp'
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
  [./stress_zz]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./strain_zz]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./peeq]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./pe11]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./pe22]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./pe33]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[AuxKernels]
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
  [./pe11]
    type = RankTwoAux
    rank_two_tensor = plastic_strain
    variable = pe11
    index_i = 0
    index_j = 0
  [../]
    [./pe22]
    type = RankTwoAux
    rank_two_tensor = plastic_strain
    variable = pe22
    index_i = 1
    index_j = 1
  [../]
  [./pe33]
    type = RankTwoAux
    rank_two_tensor = plastic_strain
    variable = pe33
    index_i = 2
    index_j = 2
  [../]
  [./eqv_plastic_strain]
    type = MaterialRealAux
    property = eqv_plastic_strain
    variable = peeq
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

  dt=0.1
  dtmax=1
  dtmin=0.1
  end_time=2.0

  nl_abs_tol = 1e-10
[]

[Postprocessors]
  [./total_stress_zz]
    type = ElementAverageValue
    variable = stress_zz
  [../]
  [./total_strain_zz]
    type = ElementAverageValue
    variable = strain_zz
  [../]
  [./Equivalent_plastic_strain]
    type = ElementAverageValue
    variable = peeq
  [../]
[]

[Outputs]
  file_base = ./tensor_mechanics_j2plasticity/out
  exodus = true
  csv = true
[]
