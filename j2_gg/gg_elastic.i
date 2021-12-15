# tar -cvf - ./test01_elastic/out_test01_elastic.e-s??[05] | pigz -9 -p 10 > out_test01_elastic.tgz


# test01_elastic, 纯弹性加载(0.8*t)，双晶box模型,界面出现扰动;
# test01_elastic_02, 纯弹性加载(5%, 50)，双晶box模型, 晶面没有出现不连续(disconnections), 界面宽度15
# test01_elastic_03, 纯弹性加载(5%, 50)，双晶box模型,  界面宽度50, 使用置换网格在求力学弱形式中,晶面没有出现不连续(disconnections),
# test01_elastic_04: 使用'if(t<40,t,40+0.02*sin(t))', ComputeIncrementalSmallStrain--ComputeFiniteStrainElasticStress, 力学-相场界面同时移动
# test01_elastic_05: 基于04,直接加入载荷4
# my_function = 'if(t<6,10*t,60+0.01*sin(t))' # 0.001s^{-1}
# my_function = '0.8*t' # 0.001s^{-1}

# my_function = 'if(t<4,t,4+0.002*sin(t))' # 0.001s^{-1}

my_filename = 'test01_elastic_07'
my_function = 'if(t<4,t,4+0.0002*sin(10*pi*t))' # 0.001s^{-1}

my_wGB = 15 #15
my_num_adaptivity = 4

my_interval = 2
my_end_time = 1e3

my_time_scale = 1.0e-9
my_length_scale = 1.0e-9
my_pressure_scale = 1.0e6

my_xmax = 3.0e2 # 3.0e3
my_ymax = 1.0e2 # 1.0e3
my_radius = 1.0e2 # 1.0e3
my_nx = 200 #200 # 50
my_ny = 50 #50 # 20

[Functions]
  [./timestep_fn]
    type = ParsedFunction
    value = 'if(t<4,0.2,5)'
  [../]
[]

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = ${my_nx}
  ny = ${my_ny}
  xmax = ${my_xmax}
  ymax = ${my_ymax}
  elem_type = QUAD4
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
  [../]
  [./disp_y]
  [../]
[]

[ICs]
  [./PolycrystalICs]
    [./BicrystalBoundingBoxIC]
      x1 = 0
      y1 = 0
      x2 = ${my_radius}
      y2 = ${my_ymax}
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
  [./elastic_stress11]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./elastic_stress22]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./elastic_stress12]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./VMstress]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./C1111]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./euler_angle]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./total_energy_density]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./grad_energy_density]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[Modules/TensorMechanics/Master]
  [./all]
    # displacements = 'disp_x disp_y'
    use_displaced_mesh = true
    strain = SMALL # FINITE
    incremental = true
  [../]
[]

[Kernels]
  [./PolycrystalKernel]
  [../]
  [./PolycrystalElasticDrivingForce]
  [../]
  # [./TensorMechanics]
  #   displacements = 'disp_x disp_y'
  #   use_displaced_mesh = true
  # [../]
[]

[AuxKernels]
  [./bnds_aux]
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
  [./elastic_stress11]
    type = RankTwoAux
    variable = elastic_stress11
    rank_two_tensor = stress
    index_i = 0
    index_j = 0
    execute_on = timestep_end
  [../]
  [./elastic_stress22]
    type = RankTwoAux
    variable = elastic_stress22
    rank_two_tensor = stress
    index_i = 1
    index_j = 1
    execute_on = timestep_end
  [../]
  [./elastic_stress12]
    type = RankTwoAux
    variable = elastic_stress12
    rank_two_tensor = stress
    index_i = 0
    index_j = 1
    execute_on = timestep_end
  [../]
  [./VMstress]
    type = RankTwoScalarAux
    rank_two_tensor = stress
    scalar_type = VonMisesStress
    variable = VMstress
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
  [./euler_angle]
    type = OutputEulerAngles
    variable = euler_angle
    euler_angle_provider = euler_angle_file
    grain_tracker = grain_tracker
    output_euler_angle = 'phi1'
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
  [./top_displacement]
    type = FunctionDirichletBC
    variable = disp_y
    boundary = top
    function = ${my_function}
  [../]
  [./x_anchor]
    type = DirichletBC
    variable = disp_x
    boundary = 'left'
    value = 0.0
  [../]
  [./y_anchor]
    type = DirichletBC
    variable = disp_y
    boundary = 'bottom'
    value = 0.0
  [../]
[]



[Materials]
  [./Copper]
    type = GBEvolution
    block = 0
    T = 500 # K
    wGB = ${my_wGB} # nm
    GBmob0 = 2.5e-6 #m^4/(Js) from Schoenfelder 1997
    Q = 0.23 #Migration energy in eV
    GBenergy = 0.708 #GB energy in J/m^2
    length_scale = ${my_length_scale}
    time_scale = ${my_length_scale}
    outputs = my_exodus
    output_properties = 'kappa_op L mu gamma_asymm sigma M_GB l_GB'
  [../]
  [./ElasticityTensor]
    type = ComputePolycrystalElasticityTensor
    grain_tracker = grain_tracker
    length_scale = ${my_length_scale}
    pressure_scale = ${my_pressure_scale}
  [../]
  [./stress]
    type = ComputeFiniteStrainElasticStress # ComputeLinearElasticStress
    block = 0
  [../]
  # [./strain]
  #   type = ComputeSmallStrain
  #   block = 0
  #   displacements = 'disp_x disp_y'
  # [../]
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

[UserObjects]
  [./euler_angle_file]
    type = EulerAngleFileReader
    file_name = test.tex
  [../]
  [./grain_tracker]
    type = GrainTrackerElasticity
    connecting_threshold = 0.05
    compute_var_to_feature_map = true
    flood_entity_type = elemental
    execute_on = 'initial timestep_begin'

    euler_angle_provider = euler_angle_file
    fill_method = symmetric9
    # C_ijkl = '1.27e5 0.708e5 0.708e5 1.27e5 0.708e5 1.27e5 0.7355e5 0.7355e5 0.7355e5'

    C_ijkl = '1.684e5 1.214e5 1.214e5 1.684e5 1.214e5 1.684e5 0.75e5 0.75e5 0.75e5'
    outputs = none
  [../]
[]

[Postprocessors]
  [./dt]
    type = TimestepSize
  [../]
  [./dofs]
    type = NumDOFs
  [../]
  [./run_time]
    type = PerfGraphData
    section_name = "Root"
    data_type = total
  [../]
  [./gr0_area]
    type = ElementIntegralVariablePostprocessor
    variable = gr0
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
  [./epsilo22_av]
    type = ElementAverageValue
    variable = elastic_strain22
  [../]
  [./sigma22_av]
    type = ElementAverageValue
    variable = elastic_stress22
  [../]
  [./VMstress_av]
    type = ElementAverageValue
    variable = VMstress
  [../]
[]

[Preconditioning]
  [./SMP]
   type = SMP
   coupled_groups = 'gr0,gr1 disp_x,disp_y'
  [../]
[]

[Executioner]
  type = Transient

  solve_type = 'PJFNK'
  petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart -pc_hypre_boomeramg_strong_threshold'
  petsc_options_value = 'hypre boomeramg 31 0.7'

  l_max_its = 30
  l_tol = 1e-4
  nl_max_its = 30
  nl_rel_tol = 1e-9

  start_time = 0.0
  end_time = ${my_end_time}

  [./TimeStepper]
    type = FunctionDT
    function = timestep_fn
    min_dt = 0.1
  [../]

  # [./TimeStepper]
  #   type = IterationAdaptiveDT
  #   dt = 1.5
  #   growth_factor = 1.2
  #   cutback_factor = 0.8
  #   optimal_iterations = 8
  # [../]

  # [./Adaptivity]
  #   initial_adaptivity = ${my_num_adaptivity}
  #   refine_fraction = 0.7
  #   coarsen_fraction = 0.1
  #   max_h_level = ${my_num_adaptivity}
  # [../]
[]

[Outputs]
  file_base = ./${my_filename}/out_${my_filename}
  csv = true
  [./my_exodus]
    type = Exodus
    interval = ${my_interval} # The interval at which time steps are output
    # sync_times = '10 50 100 500 1000 5000 10000 50000 100000'
    # sync_only = true
    sequence = true
  [../]
  [./pgraph]
    type = PerfGraphOutput
    execute_on = 'initial timestep_end final'  # Default is "final"
    level = 2                     # Default is 1
    heaviest_branch = true        # Default is false
    heaviest_sections = 7         # Default is 0
  [../]
  [./my_console]
    type = Console
    # output_linear = false
    # output_screen = false
    # interval = 1
  [../]
  execute_on = 'timestep_end'
[]
