[Mesh]
  # file = cyl-tet.e
  type = GeneratedMesh
  dim = 2
  nx = 10
  ny = 40
  xmax = 40
  ymax = 10
[]

[Variables]
  [./diffused]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[Kernels]
  [./diff]
    type = Diffusion
    variable = diffused
  [../]

  # Include our time derivative here
  [./euler]
    # type = ExampleTimeDerivative
    type = TimeDerivative
    variable = diffused
    # time_coefficient = 20.0
  [../]
[]

[BCs]
  [./bottom_diffused]
    type = DirichletBC
    variable = diffused
    boundary = left
    value = 0
  [../]

  [./top_diffused]
    type = DirichletBC
    variable = diffused
    boundary = right
    value = 1
  [../]
[]

# Transient (time-dependent) details for simulations go here:
[Executioner]
  type = Transient   # Here we use the Transient Executioner (instead of steady)
  solve_type = 'PJFNK'
  # num_steps = 75 # Run for 75 time steps, solving the system each step.
  # dt = 1 # each time step will have duration "1"
  start_time = 0.0
  end_time = 2000

  [TimeStepper]
    type = IterationAdaptiveDT
    cutback_factor = 0.9
    dt = 10.0
    growth_factor = 1.5
    optimal_iterations = 7
  []
[]

[Outputs]
  file_base = diffusion_1
  execute_on = 'timestep_end'
  exodus = true
[]
