#
# KKS simple example in the split form
#

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 150
  ny = 15
  nz = 0
  xmin = -25
  xmax = 25
  ymin = -2.5
  ymax = 2.5
  zmin = 0
  zmax = 0
  elem_type = QUAD4
[]

[AuxVariables]
  [./Fglobal]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[Variables]
  # order parameter
  [./eta]
    order = FIRST
    family = LAGRANGE
  [../]

  # hydrogen concentration
  [./c]
    order = FIRST
    family = LAGRANGE
  [../]

  # chemical potential
  [./w]
    order = FIRST
    family = LAGRANGE
  [../]

  # Liquid phase solute concentration
  [./cl]
    order = FIRST
    family = LAGRANGE
    initial_condition = 0.5
  [../]
  # Solid phase solute concentration
  [./cs]
    order = FIRST
    family = LAGRANGE
    initial_condition = 0.5
  [../]
[]

[Functions]
  [./ic_func_eta]
    type = ParsedFunction
    value = 0.5*(1.0-tanh((x)/sqrt(2.0)))
    # $0.5*(1.0-tanh(x/ \sqrt{0.2})$ 
    # 设置液相和固相$\eta$的初始条件
    # 左边界为0.5,之后急速递减到0
  [../]
  [./ic_func_c]
    type = ParsedFunction
    value = '0.5*(0.5*(1.0-tanh(x/sqrt(2.0))))^3*(6*(0.5*(1.0-tanh(x/sqrt(2.0))))^2-15*(0.5*(1.0-tanh(x/sqrt(2.0))))+10)+0.1*(1-(0.5*(1.0-tanh(x/sqrt(2.0))))^3*(6*(0.5*(1.0-tanh(x/sqrt(2.0))))^2-15*(0.5*(1.0-tanh(x/sqrt(2.0))))+10))'
    # 设置组分$c$的初始条件
    # 左边界为0.5,之后急速递减到0.1
  [../]
[]

[ICs]
  [./eta]
    variable = eta
    type = FunctionIC
    function = ic_func_eta
  [../]
  [./c]
    variable = c
    type = FunctionIC
    function = ic_func_c
  [../]
[]

# [BCs]
#   [./left_c]
#     type = DirichletBC
#     variable = 'c'
#     boundary = 'left'
#     value = 0.5
#   [../]
#   [./left_eta]
#     type = DirichletBC
#     variable = 'eta'
#     boundary = 'left'
#     value = 0.5
#   [../]
# []

[Materials]
  # Free energy of the liquid
  [./fl]
    type = DerivativeParsedMaterial
    f_name = fl
    args = 'cl'
    function = '(0.1-cl)^2'
    # 0.1--equilibrium concentrations of the liquid phase
  [../]

  # Free energy of the solid
  [./fs]
    type = DerivativeParsedMaterial
    f_name = fs
    args = 'cs'
    function = '(0.9-cs)^2'
    # 0.9--equilibrium concentrations of the liquid phase
  [../]

  # h(eta)
  [./h_eta]
    type = SwitchingFunctionMaterial
    # Helper material to provide h(η)
    # 3*eta^2-2*eta^3 Simple
    # eta^3*(6*eta^2-15*eta+10) HIGH
    h_order = HIGH
    eta = eta
  [../]

  # g(eta)
  [./g_eta]
    type = BarrierFunctionMaterial
    # Helper material to provide g(η) 
    # eta^2(1-eta)^2 simple
    # eta^2*(1-eta^2)^2 High
    g_order = SIMPLE
    eta = eta
  [../]

  # constant properties
  [./constants]
    type = GenericConstantMaterial
    prop_names  = 'M   L   eps_sq'
    prop_values = '0.7 0.7 1.0  '
  [../]
[]

[Kernels]
  # enforce c = (1-h(eta))*cl + h(eta)*cs
  [./PhaseConc]
    type = KKSPhaseConcentration
    ca       = cl
    variable = cs
    c        = c
    eta      = eta
  [../]

  # enforce pointwise equality of chemical potentials
  [./ChemPotSolute]
    type = KKSPhaseChemicalPotential
    variable = cl
    cb       = cs
    fa_name  = fl
    fb_name  = fs
  [../]

  #
  # Cahn-Hilliard Equation
  #
  [./CHBulk]
    type = KKSSplitCHCRes
    variable = c
    ca       = cl
    fa_name  = fl
    w        = w
  [../]

  [./dcdt]
    type = CoupledTimeDerivative
    variable = w
    v = c
  [../]
  [./ckernel]
    type = SplitCHWRes
    mob_name = M
    variable = w
  [../]

  #
  # Allen-Cahn Equation
  #
  [./ACBulkF]
    type = KKSACBulkF
    variable = eta
    fa_name  = fl 
    fb_name  = fs
    w        = 1.0
    args = 'cl cs'
  [../]
  [./ACBulkC]
    type = KKSACBulkC
    variable = eta
    ca       = cl
    cb       = cs
    fa_name  = fl
  [../]
  [./ACInterface]
    type = ACInterface
    variable = eta
    kappa_name = eps_sq
  [../]
  [./detadt]
    type = TimeDerivative
    variable = eta
  [../]
[]

[AuxKernels]
  [./GlobalFreeEnergy]
    variable = Fglobal
    type = KKSGlobalFreeEnergy
    fa_name = fl
    fb_name = fs
    w = 1.0
  [../]
[]

[Executioner]
  type = Transient
  solve_type = 'PJFNK'

  petsc_options_iname = '-pc_type -sub_pc_type -sub_pc_factor_shift_type'
  petsc_options_value = 'asm      ilu          nonzero'

  l_max_its = 100
  nl_max_its = 100
  nl_abs_tol = 1e-10
  start_time = 0
  end_time = 4000
  # dt = 4.0
  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 4 # Initial time step.  In this simulation it changes.
    optimal_iterations = 6 # Time step will adapt to maintain this number of nonlinear iterations
  [../]
[]

#
# Precondition using handcoded off-diagonal terms
#
[Preconditioning]
  [./full]
    type = SMP
    full = true
  [../]
[]

[Postprocessors]
  [./dofs]
    type = NumDOFs
  [../]
  [./integral]
    type = ElementL2Error
    variable = eta
    function = ic_func_eta
  [../]
[]

[VectorPostprocessors]
  [./c]
    type =  LineValueSampler
    start_point = '-25 0 0'
    end_point = '25 0 0'
    variable = c
    num_points = 151 # The number of points to sample along the line
    sort_by =  id
    execute_on = timestep_end
    # none, initial, linear, nonlinear, timestep_end, timestep_begin, final, custom
  [../]
  [./eta]
    type =  LineValueSampler
    start_point = '-25 0 0'
    end_point = '25 0 0'
    variable = eta
    num_points = 151
    sort_by =  id
    execute_on = timestep_end
  [../]
[]

[Outputs]
  [./exodus]
    type = Exodus
    execute_on = timestep_end
  [../]
  exodus = true
  console = true
  gnuplot = true
  [./csv]
    type = CSV
    execute_on = final
  [../]
[]
