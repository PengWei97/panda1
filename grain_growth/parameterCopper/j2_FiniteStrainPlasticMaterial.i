my_num_adaptivity = 3
my_filename = 'para01'

my_end_time = 1e3
my_time_scale = 1.0e-9
my_length_scale = 1.0e-7

# my_yield_0 = 700 # MPa
my_xmax = 12.0e1 # 3.0e3
my_ymax = 12.0e2 # 1.0e3
my_radius = 5.0e1 # 5 miu m 1.0e3

my_xpoint = 12.0e1 # 3.0e3
my_ypoint = 12.0e2 # 1.0e3

my_nx = 50 #200 # 50
my_ny = 50 #50 # 20

my_wGB = 0.5
my_num_adaptivity = 3

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = ${my_nx}
  ny = ${my_ny}
  xmax = ${my_xmax}
  ymax = ${my_ymax}
  elem_type = QUAD4
  uniform_refine = 2
[]

[GlobalParams]
  op_num = 2
  var_name_base = gr
[]

[Variables]
  [./PolycrystalVariables]
  [../]
[]

[ICs]
  [./PolycrystalICs]
    [./BicrystalCircleGrainIC	]
      radius = ${my_radius}
      x = ${my_xmax/2}
      y = ${my_ymax/2}
      int_width = 15 
    [../]
  [../]
[]

[AuxVariables]
  [./bnds]
    order = FIRST
    family = LAGRANGE
  [../]
  [./unique_grains]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./var_indices]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[Kernels]
  [./PolycrystalKernel]
  [../]
[]

[AuxKernels]
  [./bnds_aux]
    type = BndsCalcAux
    variable = bnds
    execute_on = timestep_end
  [../]
  [./unique_grains]
    type = FeatureFloodCountAux
    variable = unique_grains
    flood_counter = grain_tracker
    execute_on = 'initial timestep_begin'
    field_display = UNIQUE_REGION
  [../]
  [./var_indices]
    type = FeatureFloodCountAux
    variable = var_indices
    flood_counter = grain_tracker
    execute_on = 'initial timestep_begin'
    field_display = VARIABLE_COLORING
  [../]
[]

[BCs]
  [./Periodic]
    [./All]
      auto_direction = 'x'
      variable = 'gr0 gr1'
    [../]
  [../]
[]


[Materials]
  [./Copper]
    type = GBEvolution
    block = 0
    T = 500 # K
    wGB = 15 # nm
    GBmob0 = 2.5e-6 #m^4/(Js) from Schoenfelder 1997
    Q = 0.23 # Migration energy in eV
    GBenergy = 0.708 # GB energy in J/m^2
    time_scale = ${my_time_scale}
    length_scale = ${my_length_scale}
  [../]
  [./local_free_energy]
    type = DerivativeParsedMaterial
    f_name= f_chem
    args = 'gr0 gr1'
    material_property_names = 'mu gamma_asymm'
    function = 'mu*(gr0^4/4.0 - gr0^2/2.0 + gr1^4/4.0 - gr1^2/2.0 + gamma_asymm*gr0^2*gr1^2+1.0/4.0)'
    derivative_order = 2
    enable_jit = true
    outputs = my_exodus
    output_properties = 'f_chem df_chem/dgr0 df_chem/dgr1'
  [../]
[]


[Postprocessors]
  [./timestep_pp]
    type = FunctionValuePostprocessor
    function = timestep_fn
  [../]
  [./dt]
    type = TimestepSize
  [../]
  [./gr0_area]
    type = ElementIntegralVariablePostprocessor
    variable = gr0
  [../]
  [./active_time]           # Time computer spent on simulation
    type = PerfGraphData
    section_name = "Root"
    data_type = total
  [../]
  [./F_chem]
    type = ElementIntegralMaterialProperty
    mat_prop = f_chem
    # outputs = csv
  [../]
[]

[Preconditioning]
  [./SMP]
    type = SMP
    full = true
  [../]
[]

[Executioner]
  type = Transient

  scheme = bdf2 # implicit-euler # TimeIntegrator,bdf2, newmark-beta
  solve_type = PJFNK # NEWTON, PJFNK: Preconditioned Jacobian-Free Newton Krylov, Jacobian-Free Newton Krylov: Jacobian-Free Newton Krylov
  # JFNK: the Krylov solver does not require the Jacobian itself but instead the action of the Jacobian on a vector

  # solve_type = 'NEWTON'
  # petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart -pc_hypre_boomeramg_strong_threshold'
  # petsc_options_value = 'hypre boomeramg 31 0.7'

  # petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart -mat_mffd_type'
  # petsc_options_value = 'hypre    boomeramg      101                ds'

  petsc_options_iname = '-pc_type  -sub_pc_type '
  petsc_options_value = 'asm       lu' # ilu

  # petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart -pc_hypre_boomeramg_strong_threshold'
  # petsc_options_value = 'hypre boomeramg 31 0.7'
  # Values of PETSc name/value pairs

  # petsc_options = '-snes_ksp_ew'
  # petsc_options_iname = '-ksp_gmres_restart'
  # petsc_options_value = '101'

  steady_state_detection = false # teady state and stop the solve when it's reached

  # petsc_options_iname = '-pc_type' 
  # petsc_options_value = 'lu' 
  # Values of PETSc name/value pairs

  # automatic_scaling = true

  l_max_its = 30 # Max Linear Iterations
  l_tol = 1e-4 # Linear Tolerance

  nl_max_its = 30 # Max Nonlinear Iterations
  nl_rel_tol = 1e-9 # Nonlinear Relative Tolerance

  dtmin = 2.0e-6 # The minimum timestep size in an adaptive run
  # dtmax = 0.2
  start_time = 0.0
  end_time = ${my_end_time}
  # dt = 'if(t<4,0.2,1)'
  

  # [./TimeStepper]
  #   type = FunctionDT
  #   function = timestep_fn
  #   min_dt = 0.1
  # [../]


  # [./TimeStepper]
  #   type = CSVTimeSequenceStepper
  #   file_name = timesequence.csv
  #   column_name = time1
  # [../]

  # # num_steps = 30
  # dt = 0.1
  [./TimeStepper]
    type = IterationAdaptiveDT 
    # SolutionTimeAdaptiveDT: Compute simulation timestep based on actual solution time.
    dt = 0.2
    growth_factor = 1.2
    cutback_factor = 0.8
    optimal_iterations = 8
    # timestep_limiting_postprocessor = timestep_pp
    timestep_limiting_function = matl_ts_min
    time_t = '0 4 10'
    time_dt = '0.2 1 5'
  [../]
  [./Adaptivity]
    initial_adaptivity = ${my_num_adaptivity}
    refine_fraction = 0.7
    coarsen_fraction = 0.1
    max_h_level = ${my_num_adaptivity}
  [../]
[]

[Outputs]
  file_base = ./${my_filename}/out_${my_filename} 
  execute_on = 'timestep_end'
  [./my_exodus]
    type = Exodus
    interval = 2
    # append_date = true
    # append_date_format = '%Y-%m-%d'
    # discontinuous = true
    # sequence = true
  [../] 
  csv = true
  [./my_console]
    type = Console
    output_linear = false
    # output_screen = false
    # interval = 5
  [../]
  [./pgraph]
    type = PerfGraphOutput
    execute_on = 'initial timestep_end'  # Default is "final"
    level = 2                     # Default is 1
    heaviest_branch = true        # Default is false
    heaviest_sections = 7         # Default is 0
    # interval = 5
  [../]
[]

