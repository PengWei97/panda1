# sigma_y0: 2e3;
# hardFactor: 

# 'test02_j2_old_00', 使用FiniteStrainPlasticMaterial

# 'test02_j2_old_01', 使用Test1FiniteStrainPlasticMaterial(其中每个晶粒内硬化指数都是材料参数，并不与序参数关联),使用时间和网格自适应技术, hypre
# *** ERROR *** 0.497
# The following error occurred in the object "fplastic", of type "Test1FiniteStrainPlasticMaterial".
# Constitutive failure

# test02_j2_old_02,同1,使用lu
# *** ERROR *** 0.597
# The following error occurred in the object "fplastic", of type "Test1FiniteStrainPlasticMaterial".Constitutive failure

# test02_j2_old_03,同1,使用lu,使用 dtmax = 0.2 的约束
# 可以正常运算下去,不过在加载阶段计算非常慢, 计算40s的模拟时间需要花费2.67h
# 晶界进行了缓慢的迁移,动力学呈现线性变化, 但是计算消耗太大. 计算到1644s,需要花费13.6h

# test02_j2_old_04,同1,使用lu,使用 dtmax = function and scheme = BDF2 的约束
# 可以正常运算下去,不过在加载阶段计算非常慢, 计算40s的模拟时间需要花费2.56h

# test02_j2_old_05,同1,使用lu,使用 dtmax = function and scheme = BDF2 的约束
# 调试, 加了一个Preconditioning。并将sigma_y0=3e3;
# 可以正常运算下去,不过在加载阶段计算非常慢, 计算40s的模拟时间需要花费2.56h


# test02_j2_old_06,使用Test2FiniteStrainPlasticMaterial,使用lu,使用 dtmax = function and scheme = BDF2 的约束, value_0 = 3.0e3;


# test02_j2_old_07,使用Test3FiniteStrainPlasticMaterial(硬化指数梯度化处理),使用lu,使用 dtmax = function and scheme = BDF2 的约束, value_0 = 3.0e3; good
# 模拟完成，在9e3的时候晶界已经迁移到末端
# 曲线图效果可以，但是后面的应力应变分布云图不太行

# test02_j2_old_08,使用Test3FiniteStrainPlasticMaterial(硬化指数梯度化处理),使用lu,使用 dtmax = function and scheme = BDF2 的约束, value_0 = 2.0e3;没有网格自适应
# 在t=13.879539006989时，计算不收敛

# test02_j2_old_09,使用Test3FiniteStrainPlasticMaterial(硬化指数梯度化处理),使用lu,使用 dtmax = function and scheme = BDF2 的约束, value_0 = 2.0e3,打开网格自适应
# 计算非常慢，算了12.89h才算到19.56s，太慢啦！！

# test02_j2_old_10,使用使用Test2FiniteStrainPlasticMaterial(只考虑J2中的弹性能)，使用lu, 

# test02_j2_old_11: 使用ComputeIncrementalSmallStrain,strain = small, value_0 = 1.0e4, 关闭网格自适应，细化网格
# 计算到0.999超过，计算不收敛

# test02_j2_old_12: 使用ComputeIncrementalSmallStrain,strain = small, value_0 = 1.0e4, 关闭网格自适应，细化网格，缩小模拟面积
# 算到最后晶界弥散化
# 沿x方向的应变还是没有移动
  # 求解器的区别：lu 或者 hypre
  # 时间自适应和网格自适应

# test02_j2_old_13: 使用ComputeIncrementalSmallStrain,strain = small, value_0 = 1.0e4, 关闭网格自适应，细化网格，缩小模拟面积,本构修改为sigma=C_ijkl * (strain_p-delta_strain);
# 力场和相场的界面都移动了
# 沿x方向的应变还是没有移动

# test02_j2_old_14: 使用ComputeIncrementalSmallStrain,strain = small, value_0 = 1.0e4, 关闭网格自适应，细化网格，缩小模拟面积,本构修改为sigma=C_ijkl * (strain_p-delta_strain), 并打开returnMap部分；
# 可以计算，没有问题

# test02_j2_old_15: 使用ComputeIncrementalSmallStrain,strain = small, value_0 = 1.0e4, 关闭网格自适应，细化网格，缩小模拟面积,本构修改为sigma=C_ijkl * (strain_p-delta_strain), 并打开returnMap部分；
# 力场和相场的界面都移动了，界面移动太慢

# test02_j2_old_16: 使用ComputeIncrementalSmallStrain,strain = small, value_0 = 1.0e4, 关闭网格自适应，细化网格，缩小模拟面积,本构修改为sigma=C_ijkl * (strain_p-delta_strain), 并打开returnMap部分；
# 力场和相场的界面都移动了,但是力场界面后面没有跟着移动


# test02_j2_old_16: 使用ComputeIncrementalSmallStrain,strain = small, value_0 =2.0e4, 关闭网格自适应，细化网格，缩小模拟面积,本构修改为sigma=C_ijkl * (strain_total+delta_strain-strain_plastic), 并打开returnMap部分,计算域在弹性阶段；
# 力场和相场的界面都移动了,但是力场界面后面没有跟着移动

# test02_j2_old_17: 使用ComputeIncrementalSmallStrain,strain = small, value_0 = 3.0e3, 关闭网格自适应，细化网格，缩小模拟面积,本构修改为sigma=C_ijkl * (strain_total+delta_strain-strain_plastic), 并打开returnMap部分,计算域在弹性阶段；
# 力场和相场的界面都移动了,但是力场界面后面没有跟着移动

# test02_j2_old_18: 使用ComputeIncrementalSmallStrain,strain = small, value_0 = 3.0e3, 关闭网格自适应，细化网格，缩小模拟面积,本构修改为sigma=C_ijkl * (strain_total+delta_strain-strain_plastic), 并打开returnMap部分,计算域在弹性阶段；
# 力场和相场的界面都移动了,但是力场界面后面没有跟着移动
# 需要进一步分析

# test02_j2_old_19: 使用ComputeIncrementalSmallStrain,strain = small, value_0 = 3.0e3, 关闭网格自适应，细化网格，缩小模拟面积,本构修改为sigma=C_ijkl * (strain_total+delta_strain-strain_plastic), 并打开returnMap部分,计算域在弹性阶段；
# 力场和相场的界面都移动了,但是力场界面后面没有跟着移动
# 需要进一步分析

my_filename = 'test02_j2_old_16'
# my_function = 'if(t<40,t,40+0.02*sin(t))' # 0.001s^{-1}

my_function = 'if(t<4,t,4+0.0002*sin(10*pi*t))' # 0.001s^{-1}
my_end_time = 1e3
# my_yield_0 = 700 # MPa
my_xmax = 3.0e2 # 3.0e3
my_ymax = 1.0e2 # 1.0e3
my_radius = 1.0e2 # 1.0e3
my_nx = 200 #200 # 50
my_ny = 50 #50 # 20

my_time_scale = 1.0e-9
my_length_scale = 1.0e-9
my_pressure_scale = 1.0e6

my_yield_strength_init = 2e3

[Functions]
  [./timestep_fn]
    type = ParsedFunction
    value = 'if(t<4,0.2,5)'
  [../]
  [./dts]
    type = PiecewiseLinear
    x = '0 4 1000'
    y = '0 0.2 4'
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
  # uniform_refine = 2
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
  [./strain_xx]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./strain_xy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./strain_yy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./elastic_strain11]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./elastic_strain12]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./elastic_strain22]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./plastic_strain11]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./plastic_strain12]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./plastic_strain22]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_xx]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_xy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_yy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./vm_stress]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./eq_pl_strain]
    order = CONSTANT
    family = MONOMIAL
  [../]

  # [./p_internal_parameter]
  #   order = CONSTANT
  #   family = MONOMIAL
  # [../]
  [./HardFactor]
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
  [./C1111]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./euler_angle]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[Kernels]
  [./PolycrystalKernel]
  [../]
  [./PolycrystalElasticDrivingForce]
    # ACGrGrElasticDrivingForce
    # GGPolycrystalElasticDrivingForce
    op_num = 2
    var_name_base = gr
      # GGACGrGrElasticDrivingForce
  [../]
  # [./ACGrGrPlasticDrivingForce1]
  #   type = ACGrGrPlasticDrivingForce
  #   D_hard_factor_name = dhard_factor/dgr0
  #   variable = gr0
  # [../]
  # [./ACGrGrPlasticDrivingForce2]
  #   type = ACGrGrPlasticDrivingForce
  #   D_hard_factor_name = dhard_factor/dgr1
  #   variable = gr1
  # [../]
[]

[Modules/TensorMechanics/Master]
  [./all]
    # strain = FINITE
    # displacements = 'disp_x disp_y'  
    # use_displaced_mesh = true
    # add_variables = true
    # use_automatic_differentiation = true
    strain = SMALL # FINITE
    incremental = true
  [../]
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
  [./elastic_strain12]
    type = RankTwoAux
    variable = elastic_strain11
    rank_two_tensor = elastic_strain
    index_i = 0
    index_j = 1
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
  [./plastic_strain11]
    type = RankTwoAux
    variable = plastic_strain11
    rank_two_tensor = plastic_strain
    index_i = 0
    index_j = 0
    execute_on = timestep_end
  [../]
  [./plastic_strain12]
    type = RankTwoAux
    variable = plastic_strain12
    rank_two_tensor = plastic_strain
    index_i = 0
    index_j = 1
    execute_on = timestep_end
  [../]
  [./plastic_strain22]
    type = RankTwoAux
    variable = plastic_strain22
    rank_two_tensor = plastic_strain
    index_i = 1
    index_j = 1
    execute_on = timestep_end
  [../]
  [./strain_xx]
    type = RankTwoAux
    rank_two_tensor = total_strain
    variable = strain_xx
    index_i = 0
    index_j = 0
  [../]
  [./strain_xy]
    type = RankTwoAux
    rank_two_tensor = total_strain
    variable = strain_xy
    index_i = 0
    index_j = 1
  [../]
  [./strain_yy]
    type = RankTwoAux
    rank_two_tensor = total_strain
    variable = strain_yy
    index_i = 1
    index_j = 1
  [../]
  [./stress_xx]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_xx
    index_i = 0
    index_j = 0
  [../]
  [./stress_xy]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_yy
    index_i = 0
    index_j = 1
  [../]
  [./stress_yy]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_yy
    index_i = 1
    index_j = 1
  [../]
  # [./intnl]
  #   type = MaterialStdVectorAux
  #   index = 0
  #   property = plastic_internal_parameter
  #   variable = intnl
  # [../]
  # [./eq_pl_strain]
  #   type = RankTwoScalarAux
  #   rank_two_tensor = plastic_strain
  #   scalar_type = EffectiveStrain
  #   variable = eq_pl_strain
  # [../]
  [./vm_stress]
    type = RankTwoScalarAux
    rank_two_tensor = stress
    scalar_type = VonMisesStress
    variable = vm_stress
  [../]
  [./eq_pl_strain]
    type = MaterialRealAux
    # index = 0
    property = eqv_plastic_strain
    variable = eq_pl_strain
  [../]
  [./HardFactor]
    type = MaterialRealAux
    variable = HardFactor
    property = hard_factor  
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
[]

[BCs]
  [./top_displacement]
    type = FunctionDirichletBC
    variable = disp_y
    boundary = top
    function = ${my_function}
    # use_displaced_mesh = true
    # function = 50
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
    wGB = 15 # nm
    GBmob0 = 2.5e-6 #m^4/(Js) from Schoenfelder 1997
    Q = 0.23 # Migration energy in eV
    GBenergy = 0.708 # GB energy in J/m^2
    time_scale = ${my_time_scale}
    length_scale = ${my_length_scale}
  [../]
  [./ElasticityTensor]
    type = ComputePolycrystalElasticityTensor
    grain_tracker = grain_tracker
  [../]
  [./fplastic]
    type = Test2FiniteStrainPlasticMaterial #设置屈服函数
    # 选择: ComputeMultipleInelasticStress
    # implements rate-independent associative J2 plasticity 
    # with isotropic hardening in the finite-strain framework.
    block = 0
    grain_tracker = grain_tracker
    yield_stress='0. 700. 0.05 700. 0.1 700. 0.38 700. 0.95 700. 2. 700.'
    # outputs = my_exodus
    # output_properties = 'dhard_factor/dgr0 dhard_factor/dgr1'
    # length_scale = ${my_length_scale}
    # pressure_scale = ${my_pressure_scale}
    # yield_strength_init = ${my_yield_strength_init}
  [../]
  # [./strain]
  #   type = ComputeFiniteStrain
  #   block = 0
  #   displacements = 'disp_x disp_y'
  # [../]
  [./elastic_free_energy]
    type = ElasticEnergyMaterial
    f_name = f_elastic
    block = 0
    args = 'gr0 gr1'
    outputs = my_exodus
    output_properties = 'f_elastic df_elastic/dgr0 df_elastic/dgr1'
  [../]
  # [./plastic_free_energy]
  #   type = PlasticEnergyMaterial
  #   f_name = f_plastic
  #   block = 0
  #   args = 'gr0 gr1'
  #   outputs = my_exodus
  #   output_properties = 'f_plastic df_plastic/dgr0 df_plastic/dgr1'
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
    C_ijkl = '1.684e5 1.214e5 1.214e5 1.684e5 1.214e5 1.684e5 0.75e5 0.75e5 0.75e5'

    outputs = none
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
  [./epsilon_yy]
    type = PointValue
    point = '0 0 0'
    variable = strain_yy
  [../]
  [./epsilon_yy_av]
    type = ElementAverageValue
    variable = strain_yy
  [../]
  [./sigma_yy]
    type = PointValue
    point = '0 0 0'
    variable = stress_yy
  [../]
  [./sigma_yy_av]
    type = ElementAverageValue
    variable = stress_yy
  [../]
  [./VMstress]
    type = PointValue
    point = '0 0 0'
    variable = vm_stress
  [../]
  [./VMstress_av]
    type = ElementAverageValue
    variable = vm_stress
  [../]
  [./eq_pl_strain]
    type = PointValue
    point = '0 0 0'
    variable = eq_pl_strain
  [../]
  [./eq_pl_strain_av]
    type = ElementAverageValue
    variable = eq_pl_strain
  [../]
  [./HardFactor]
    type = PointValue
    point = '0 0 0'
    variable = HardFactor
  [../]
  [./HardFactor_av]
    type = ElementAverageValue
    variable = HardFactor
  [../]
  [./active_time]           # Time computer spent on simulation
    type = PerfGraphData
    section_name = "Root"
    data_type = total
  [../]
  [./F_elastic]
    type = ElementIntegralMaterialProperty
    mat_prop = f_elastic
  [../]
  # [./F_plastic]
  #   type = ElementIntegralMaterialProperty
  #   mat_prop = f_plastic
  # [../]
  [./F_chem]
    type = ElementIntegralMaterialProperty
    mat_prop = f_chem
    # outputs = csv
  [../]
[]

[Preconditioning]
  [./SMP]
   type = SMP
   coupled_groups = 'gr0,gr1 disp_x,disp_y'
  [../]
[]

#
# Precondition using handcoded off-diagonal terms
#
# [Preconditioning]
#   [./full]
#     type = SMP
#     full = true
#   [../]
# []

# [Preconditioning]
#   [./SMP]
#     type = SMP
#     off_diag_row = 'disp_x disp_y'
#     off_diag_column = 'disp_y disp_x'
#   [../]
# []

# [Preconditioning]
#   [./SMP]
#    type = SMP
#    coupled_groups = 'disp_x,disp_y'  # poly_grain_growth_2D_eldrforce.i
#   [../]
# []

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
  

  [./TimeStepper]
    type = FunctionDT
    function = timestep_fn
    min_dt = 0.1
  [../]


  # [./TimeStepper]
  #   type = CSVTimeSequenceStepper
  #   file_name = timesequence.csv
  #   column_name = time1
  # [../]

  # # num_steps = 30
  # dt = 0.1
  # [./TimeStepper]
  #   type = IterationAdaptiveDT 
  #   # SolutionTimeAdaptiveDT: Compute simulation timestep based on actual solution time.
  #   dt = 0.2
  #   growth_factor = 1.2
  #   cutback_factor = 0.8
  #   optimal_iterations = 8
  #   # timestep_limiting_postprocessor = timestep_pp
  #   timestep_limiting_function = matl_ts_min
  #   time_t = '0 4 10'
  #   time_dt = '0.2 1 5'
  # [../]
  # [./Adaptivity]
  #   initial_adaptivity = 5
  #   refine_fraction = 0.7
  #   coarsen_fraction = 0.1
  #   max_h_level = 3
  # [../]
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
    sequence = true
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

# j2_hard1_mod_optimised.i
# PowerLawHardening.i

