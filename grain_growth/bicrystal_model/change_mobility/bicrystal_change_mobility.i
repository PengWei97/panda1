my_filename = 'mGB_0012' 
my_wGB = 30
my_mobilityGB = 0.12e-13

my_xmax = 400e3
my_point = 200e3
my_radius = 180e3

my_length_scale = 1.0e-9
my_time_scale = 1.0e-3

my_num_adaptivity = 3
my_interval = 2
my_loading = 0
my_uniform_refine = 2

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 40
  ny = 40
  xmax = ${my_xmax}
  ymax = ${my_xmax}
  elem_type = QUAD4
  uniform_refine = ${my_uniform_refine}
[]

[GlobalParams]
  op_num = 2
  var_name_base = gr
  displacements = 'disp_x disp_y'
[]

[Variables]
  [./PolycrystalVariables]
  [../]
  [./disp_x]
    order = FIRST
    family = LAGRANGE
  [../]
  [./disp_y]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[UserObjects]
  [./euler_angle_file]
    type = EulerAngleFileReader
    file_name = grn_2_rand_2D.tex
  [../]
  [./grain_tracker]
    type = GrainTrackerElasticity
    threshold = 0.2
    compute_var_to_feature_map = true
    execute_on = 'initial timestep_begin'
    flood_entity_type = ELEMENTAL

    C_ijkl = '1.27e5 0.708e5 0.708e5 1.27e5 0.708e5 1.27e5 0.7355e5 0.7355e5 0.7355e5' # FCC Cu
    # C_ijkl = '2.31e5 1.347e5 1.347e5 2.31e5 1.347e5 2.31e5 1.164e5 1.164e5 1.164e5'    # BCC Fe
    # C_ijkl = '1.94e5 0.655e5 0.698e5 1.94e5 0.698e5 1.98e5 0.4627e5 0.4627e5 0.6435e5' # Titanium,2,0Pa，可行
    # C_ijkl = '194.305e3	65.597e3	69.870e3	194.305e3 69.870e3 198.907e3 	46.270e3	46.270e3 64.354e3' # Titanium,0MPa，可行
    # C_ijkl = '1.60e5 0.90e5 0.66e5 1.60e5 0.66e5 1.81e5 0.465e5 0.465e5 0.350e5' # Titanium, table
    # C_ijkl = '1.60e5 0.90e5 0.66e5 1.60e5 0.66e5 1.81e5 0.465e5 0.465e5 0.35e5' 
     
    fill_method = symmetric9
    euler_angle_provider = euler_angle_file
  [../]
[]

[ICs]
  [./PolycrystalICs]
    [./BicrystalCircleGrainIC]
      radius = ${my_radius}
      x = ${my_point}
      y = ${my_point}
      int_width = ${my_wGB}
    [../]
  [../]
[]

[AuxVariables]
  [./bnds]
    order = FIRST
    family = LAGRANGE
  [../]
  [./elastic_strain11]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./elastic_strain22]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./elastic_strain12]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress11]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress12]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress22]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./unique_grains]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./var_indices]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./vonmises_stress]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./C1111]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./euler_angle_sigma]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./euler_angle_phi]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./total_energy_density]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

# [Modules/TensorMechanics/Master]
#   [./all]
#     # displacements = 'disp_x disp_y'
#     use_displaced_mesh = true
#     strain = SMALL # FINITE
#     # incremental = true
#   [../]
# []

[Kernels]
  [./PolycrystalKernel]
  [../]
  [./PolycrystalElasticDrivingForce]
  [../]
  [./TensorMechanics]
    use_displaced_mesh = true
    displacements = 'disp_x disp_y'
  [../]
[]

[AuxKernels]
  [./BndsCalc]
    type = BndsCalcAux
    variable = bnds
    execute_on = timestep_end
  [../]
  [./elastic_strain11]
    type = RankTwoAux
    variable = elastic_strain11
    rank_two_tensor = elastic_strain
    index_i = 0
    index_j = 0
    execute_on = timestep_end
  [../]
  [./elastic_strain22]
    type = RankTwoAux
    variable = elastic_strain22
    rank_two_tensor = elastic_strain
    index_i = 1
    index_j = 1
    execute_on = timestep_end
  [../]
  [./elastic_strain12]
    type = RankTwoAux
    variable = elastic_strain12
    rank_two_tensor = elastic_strain
    index_i = 0
    index_j = 1
    execute_on = timestep_end
  [../]
  [./stress11]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress11
    index_i = 0
    index_j = 0
  [../]
  [./stress12]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress12
    index_i = 0
    index_j = 1
  [../]
  [./stress22]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress22
    index_i = 1
    index_j = 1
  [../]
  [./unique_grains]
    type = FeatureFloodCountAux
    variable = unique_grains
    execute_on = timestep_end
    flood_counter = grain_tracker
    field_display = UNIQUE_REGION
  [../]
  [./var_indices]
    type = FeatureFloodCountAux
    variable = var_indices
    execute_on = timestep_end
    flood_counter = grain_tracker
    field_display = VARIABLE_COLORING
  [../]
  [./C1111]
    type = RankFourAux
    variable = C1111
    rank_four_tensor = elasticity_tensor
    index_l = 0
    index_j = 0
    index_k = 0
    index_i = 0
    execute_on = timestep_end
  [../]
  [./vonmises_stress]
    type = RankTwoScalarAux
    variable = vonmises_stress
    rank_two_tensor = stress
    scalar_type = VonMisesStress
    execute_on = timestep_end
  [../]
  [./euler_angle_sigma]
    type = OutputEulerAngles
    variable = euler_angle_sigma
    euler_angle_provider = euler_angle_file
    grain_tracker = grain_tracker
    output_euler_angle = 'phi1'
    #  phi1, Phi, phi2
    execute_on = 'initial timestep_end'
  [../]
  [./euler_angle_phi]
    type = OutputEulerAngles
    variable = euler_angle_phi
    euler_angle_provider = euler_angle_file
    grain_tracker = grain_tracker
    output_euler_angle = 'Phi'
    #  phi1, Phi, phi2
    execute_on = 'initial timestep_end'
  [../]
  [./local_free_energy]
    type = TotalFreeEnergy
    f_name = f_chem
    variable = total_energy_density
    kappa_names = 'kappa_op kappa_op'
    interfacial_vars = 'gr0 gr1'
  [../]
[]

[BCs]
  [./Periodic]
    [./All]
      auto_direction = 'x'
      variable = 'gr0 gr1'
    [../]
  [../]
  [./top_displacement]
    type = DirichletBC
    variable = disp_y
    boundary = top
    value = ${my_loading}
    # value = 0
  [../]
  [./x_anchor]
    type = DirichletBC
    variable = disp_x
    boundary = 'left right'
    value = 0.0
  [../]
  [./y_anchor]
    type = DirichletBC
    variable = disp_y
    boundary = bottom
    value = 0.0
  [../]
[]

[Materials]
  [./Copper]
    type = GBEvolution
    block = 0
    T = 375.15 # K
    wGB = ${my_wGB} # nm
    # GBmob0 = 2.5e-6 # m^4/(Js) from Schoenfelder 1997
    GBMobility = ${my_mobilityGB} # m^4/(J · s)
    Q = 0.23 # Migration energy in eV
    GBenergy = 0.60 # GB energy in J/m^2
    length_scale = ${my_length_scale} # μm
    time_scale = ${my_time_scale} # s
  [../]
  [./ElasticityTensor]
    type = ComputePolycrystalElasticityTensor
    grain_tracker = grain_tracker
  [../]
  [./strain]
    type = ComputeSmallStrain
    block = 0
    displacements = 'disp_x disp_y'
  [../]
  [./stress]
    type = ComputeLinearElasticStress
    block = 0
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
  [./elastic_free_energy]
    type = ElasticEnergyMaterial
    f_name = f_elastic
    block = 0
    args = 'gr0 gr1'
    outputs = my_exodus
    output_properties = 'f_elastic df_elastic/dgr0 df_elastic/dgr1'
    # f_elastic MPa
    # df_elastic/dgr0--eV/nm^2
    # MPa*(length_scale^3)*pressure_scale;
  [../]
[]

[Postprocessors]
  [./dt]
    type = TimestepSize
  [../]
  [./dofs]
    type = NumDOFs
  [../]
  [./F_elastic]
    type = ElementIntegralMaterialProperty
    mat_prop = f_elastic
  [../]
  [./F_chem]
    type = ElementIntegralMaterialProperty
    mat_prop = f_chem
    # outputs = csv
  [../]
  [./run_time]
    type = PerfGraphData
    section_name = "Root"
    data_type = total
  [../]
  [./gr0area]
    type = ElementIntegralVariablePostprocessor
    variable = gr0
  [../]
  [./VMstress_av]
    type = ElementAverageValue
    variable = vonmises_stress
  [../]
[]

[Preconditioning]
  [./SMP]
    type = SMP
    coupled_groups = 'gr0,gr1,disp_x,disp_y'
  [../]
[]

[Executioner]
  type = Transient
  scheme = bdf2
  solve_type = PJFNK
  petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart -pc_hypre_boomeramg_strong_threshold'
  petsc_options_value = 'hypre boomeramg 31 0.7'
  l_tol = 1.0e-4
  l_max_its = 30
  nl_max_its = 25
  nl_rel_tol = 1.0e-7

  start_time = 0.0
  end_time = 7200

  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 1.5
    growth_factor = 1.2
    cutback_factor = 0.8
    optimal_iterations = 8
  [../]
  [./Adaptivity]
    initial_adaptivity = ${my_num_adaptivity}
    refine_fraction = 0.8 # The fraction of elements or error to refine. Should be between 0 and 1.
    coarsen_fraction = 0.05 # 0.05 Fraction of low error that will coarsened
    max_h_level = ${my_num_adaptivity}
  [../]
[]

[Outputs]
  file_base = ./${my_filename}/out_${my_filename}
  [./my_exodus]
    type = Exodus
    interval = ${my_interval} # The interval at which time steps are output
    # sync_times = '10 50 100 500 1000 5000 10000 50000 100000'
    # sync_only = true
    sequence = true
  [../]
  csv = true
  [pgraph]
    type = PerfGraphOutput
    execute_on = 'initial final'  # Default is "final"
    level = 2                     # Default is 1
    heaviest_branch = true        # Default is false
    heaviest_sections = 7         # Default is 0
  []
[]