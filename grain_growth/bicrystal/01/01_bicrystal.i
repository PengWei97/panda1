[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 20
  ny = 6
  xmax = 1000
  ymax = 1000
  elem_type = QUAD4
  uniform_refine = 1
[]

[GlobalParams]
  op_num = 2
  var_name_base = gr
  length_scale = 1.0e-9
  time_scale = 1.0e-6
  pressure_scale = 1.0e6
  # use_displaced_mesh = False
[]

[Variables]
  [./PolycrystalVariables]
     # PhaseFieldApp
     # Variables/PolycrystalVariables --> PolycrystalVariablesAction
  [../]
  [./disp_x]
  [../]
  [./disp_y]
  [../]
[]

[ICs]
  [./PolycrystalICs]
    [./BicrystalBoundingBoxIC]
     # PhaseFieldApp
     # ICs/PolycrystalICs/BicrystalBoundingBoxIC --> BicrystalBoundingBoxICAction
      x1 = 0
      y1 = 0
      x2 = 500
      y2 = 1000
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
  [./active_bounds_elemental]
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
     # registerSyntax("PolycrystalKernelAction", "Kernels/PolycrystalKernel");
     # input:_op_num = int,_var_name_base = gr
     # kernel:
          # ACGrGrPoly: Interface Allen-Cahn Kernel,
            # \mu*(\eta_1^3-\eta_1+2*\gamma*\eta_1*\eta_2^2) for bicrystal
          # ACInterface: the Gradient energy Allen-Cahn Kernel
          # TimeDerivative:
          # ACGBPoly: 保守场序参量和非保守场序参量之间的耦合项，用于表示沉淀
  [../]
  # [./PolycrystalElasticDrivingForce]
  #   # registerSyntax("PolycrystalElasticDrivingForceAction", "Kernels/PolycrystalElasticDrivingForce");
  #   # ACGrGrElasticDrivingForce
  # [../]
  [./AC_ElasticDrivingForce_gr0]
    type = ACGrGrElasticDrivingForce
    variable = gr0
    D_tensor_name = delasticity_tensor/dgr0
  [../]
  [./AC_ElasticDrivingForce_gr1]
    type = ACGrGrElasticDrivingForce
    variable = gr1
    D_tensor_name = delasticity_tensor/dgr1
  [../]
  [./TensorMechanics]
    # registerSyntax("LegacyTensorMechanicsAction", "Kernels/TensorMechanics");
    # 这个遗留操作很快就会被弃用，取而代之的是更具包容性的TensorMechanics/MasterAction。
    # LegacyTensorMechanicsAction 是一个方便的对象，它简化了张量力学系统设置的一部分。 它添加了 StressDivergence Kernels（对于当前坐标系）。
    displacements = 'disp_x disp_y'
  [../]
[]

[AuxKernels]
  [./bnds_aux]
    type = BndsCalcAux
    # 用于可视化晶界
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
  [./active_bounds_elemental]
    type = FeatureFloodCountAux
    variable = active_bounds_elemental
    field_display = ACTIVE_BOUNDS
    execute_on = 'initial timestep_begin'
    flood_counter = grain_tracker
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
    type = DirichletBC
    variable = disp_y
    boundary = top
    value = -10.0
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
    T = 500 # K
    wGB = 75 # nm
    GBmob0 = 2.5e-6 #m^4/(Js) from Schoenfelder 1997
    Q = 0.23 #Migration energy in eV
    GBenergy = 0.708 #GB energy in J/m^2
    # time_scale = 1.0e-6 # \mu s
    # length_scale = 1.0e-9 # nm
  [../]
  [./ElasticityTensor]
    type = ComputePolycrystalElasticityTensor
    grain_tracker = grain_tracker
    # length_scale = 1.0e-9
    # pressure_scale = 1.0e6
    outputs = Exodus
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
    C_ijkl = '1.27e5 0.708e5 0.708e5 1.27e5 0.708e5 1.27e5 0.7355e5 0.7355e5 0.7355e5'

    # outputs = Exodus
    # 序参数重映射算法，用于减少序参数的数目
  [../]
[]

[Postprocessors]
  [./dt]
    type = TimestepSize
  [../]
  [./gr0_area]
    type = ElementIntegralVariablePostprocessor
    variable = gr0
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
  num_steps = 30
  dt = 0.2

  [./Adaptivity]
   initial_adaptivity = 2
    refine_fraction = 0.7
    coarsen_fraction = 0.1
    max_h_level = 2
  [../]
[]

[Outputs]
  execute_on = 'timestep_end'
  [./Exodus]
    type = Exodus
    # interval = 1
    # end_step = 5
  [../]
  [./pgraph]
    type = PerfGraphOutput
    execute_on = 'initial final'  # Default is "final"
    # execute_on = 'timestep_end'
    level = 2                     # Default is 1
    heaviest_branch = true        # Default is false
    heaviest_sections = 7         # Default is 0
  []
  [./csv]
    type = CSV
    # interval = 1
  [../]
[]
