[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 25
  ny = 25
  xmax = 50
  ymax = 50
  elem_type = QUAD4
[]

[Variables]
  [./c]
    order = THIRD
    family = HERMITE
  [../]
[]

[ICs]
  [./c]
    type = SmoothSuperellipsoidIC
    variable = c
    x1 = 25.0
    y1 = 25.0
    a = 8.0
    b = 12.0
    n = 3.5
    invalue = 1.0
    outvalue = -0.8
    int_width = 4.0
  [../]
[]

[Kernels]
  [./ie_c]
    type = TimeDerivative
    variable = c
  [../]
  [./CHSolid]
    type = CHMath
    # f_loc = 1/4*(1-c^2)^2, grad(df/dc) = 3*c^2*grad_c-grad_c
    variable = c
    mob_name = M
  [../]
  [./CHInterface]
    type = CHInterface
    variable = c
    kappa_name = kappa_c
    mob_name = M
  [../]
[]

[BCs]
  [./Periodic]
    [./all]
      auto_direction = 'x y'
    [../]
  [../]
[]

[Materials]
  [./constant]
    type = GenericConstantMaterial
    prop_names  = 'M kappa_c'
    prop_values = '1.0 1.0'
  [../]
[]

[Executioner]
  type = Transient
  scheme = 'bdf2'
  solve_type = 'NEWTON'

  petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart'
  petsc_options_value = 'hypre boomeramg 31'

  l_max_its = 20
  l_tol = 1.0e-4
  nl_max_its = 40
  nl_rel_tol = 1e-9

  start_time = 0.0
  num_steps = 120
  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 2.0
    growth_factor = 1.2
    cutback_factor = 0.8
    optimal_iterations = 8
  [../]
[]

[Outputs]
  exodus = false
  [./out]
    type = Exodus
    refinements = 2
  [../]
[]
