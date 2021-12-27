my_filename = 'FCC_Cu_poly200_3D' 

[Mesh]
  type = GeneratedMesh
  dim = 3
  nx = 40
  ny = 40
  nz = 1
  xmax = 2000
  ymax = 2000
  zmax = 25
  elem_type = HEX8 # QUAD4 # TRI3
  uniform_refine = 2
[]

[GlobalParams]
  op_num = 8
  var_name_base = gr
  grain_num = 200
  displacements = 'disp_x disp_y disp_z'
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
  [./disp_z]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[UserObjects]
  [./euler_angle_file]
    type = EulerAngleFileReader
    file_name = grn_200_rand_2D.tex
  [../]
  [./voronoi]
    type = PolycrystalVoronoi
    rand_seed = 300
    coloring_algorithm = jp
    columnar_3D = true
  [../]
  [./grain_tracker]
    type = GrainTrackerElasticity
    threshold = 0.2
    compute_var_to_feature_map = true
    execute_on = 'initial timestep_begin'
    flood_entity_type = ELEMENTAL

    C_ijkl = '1.27e5 0.708e5 0.708e5 1.27e5 0.708e5 1.27e5 0.7355e5 0.7355e5 0.7355e5' # begin
    # C_ijkl = '1.684e5 1.214e5 1.214e5 1.684e5 1.214e5 1.684e5 0.754e5 0.754e5 0.754e5' # copper
    # C_ijkl = '1.60e5 0.24e5 0.06e5 1.60e5 0.06e5 1.81e5 0.465e5 0.90e5 0.68e5' # Titanium
    # C_ijkl = '1111 1122 1133 2222 2233 3333 2323 1313 1212'


    # 需要进一步理解，弹性刚度矩阵，为了各向异性，以及查找每个晶粒的欧拉角演化情况，MPa
    fill_method = symmetric9
    euler_angle_provider = euler_angle_file
  [../]
[]

[ICs]
  [./PolycrystalICs]
    [./PolycrystalColoringIC]
      polycrystal_ic_uo = voronoi
    [../]
  [../]
[]

[AuxVariables]
  [./bnds]
    order = FIRST
    family = LAGRANGE
  [../]
  # [./elastic_strain11]
  #   order = CONSTANT
  #   family = MONOMIAL
  # [../]
  # [./elastic_strain22]
  #   order = CONSTANT
  #   family = MONOMIAL
  # [../]
  # [./elastic_strain12]
  #   order = CONSTANT
  #   family = MONOMIAL
  # [../]
  # [./stress11]
  #   order = CONSTANT
  #   family = MONOMIAL
  # [../]
  # [./stress12]
  #   order = CONSTANT
  #   family = MONOMIAL
  # [../]
  # [./stress22]
  #   order = CONSTANT
  #   family = MONOMIAL
  # [../]
  # [./unique_grains]
  #   order = CONSTANT
  #   family = MONOMIAL
  # [../]
  # [./var_indices]
  #   order = CONSTANT
  #   family = MONOMIAL
  # [../]
  # [./vonmises_stress]
  #   order = CONSTANT
  #   family = MONOMIAL
  # [../]
  # [./C1111]
  #   order = CONSTANT
  #   family = MONOMIAL
  # [../]
  [./euler_angle]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[Kernels]
  [./PolycrystalKernel]
  [../]
  [./PolycrystalElasticDrivingForce]
  [../]
  [./TensorMechanics]
    use_displaced_mesh = true
    displacements = 'disp_x disp_y disp_z'
  [../]
[]

[AuxKernels]
  [./BndsCalc]
    type = BndsCalcAux
    variable = bnds
    execute_on = timestep_end
  [../]
  # [./elastic_strain11]
  #   type = RankTwoAux
  #   variable = elastic_strain11
  #   rank_two_tensor = elastic_strain
  #   index_i = 0
  #   index_j = 0
  #   execute_on = timestep_end
  # [../]
  # [./elastic_strain22]
  #   type = RankTwoAux
  #   variable = elastic_strain22
  #   rank_two_tensor = elastic_strain
  #   index_i = 1
  #   index_j = 1
  #   execute_on = timestep_end
  # [../]
  # [./elastic_strain12]
  #   type = RankTwoAux
  #   variable = elastic_strain12
  #   rank_two_tensor = elastic_strain
  #   index_i = 0
  #   index_j = 1
  #   execute_on = timestep_end
  # [../]
  # [./stress11]
  #   type = RankTwoAux
  #   rank_two_tensor = stress
  #   variable = stress11
  #   index_i = 0
  #   index_j = 0
  # [../]
  # [./stress12]
  #   type = RankTwoAux
  #   rank_two_tensor = stress
  #   variable = stress12
  #   index_i = 0
  #   index_j = 1
  # [../]
  # [./stress22]
  #   type = RankTwoAux
  #   rank_two_tensor = stress
  #   variable = stress22
  #   index_i = 1
  #   index_j = 1
  # [../]
  # [./unique_grains]
  #   type = FeatureFloodCountAux
  #   variable = unique_grains
  #   execute_on = timestep_end
  #   flood_counter = grain_tracker
  #   field_display = UNIQUE_REGION
  # [../]
  # [./var_indices]
  #   type = FeatureFloodCountAux
  #   variable = var_indices
  #   execute_on = timestep_end
  #   flood_counter = grain_tracker
  #   field_display = VARIABLE_COLORING
  # [../]
  # [./C1111]
  #   type = RankFourAux
  #   variable = C1111
  #   rank_four_tensor = elasticity_tensor
  #   index_l = 0
  #   index_j = 0
  #   index_k = 0
  #   index_i = 0
  #   execute_on = timestep_end
  # [../]
  # [./vonmises_stress]
  #   type = RankTwoScalarAux
  #   variable = vonmises_stress
  #   rank_two_tensor = stress
  #   scalar_type = VonMisesStress
  #   execute_on = timestep_end
  # [../]
  [./euler_angle]
    type = OutputEulerAngles
    variable = euler_angle
    euler_angle_provider = euler_angle_file
    grain_tracker = grain_tracker
    output_euler_angle = 'phi1'
    execute_on = 'initial timestep_end'
  [../]
[]

[BCs]
  [./Periodic]
    [./All]
      auto_direction = 'x'
      variable = 'gr0 gr1 gr2 gr3 gr4 gr5 gr6 gr7'
    [../]
  [../]
  [./top_displacement]
    type = DirichletBC
    variable = disp_y
    boundary = top
    value = 40.0
  [../]
  [./left_right]
    type = DirichletBC
    variable = disp_x
    boundary = 'left right'
    value = 0.0
  [../]
  [./bottom_y]
    type = DirichletBC
    variable = disp_y
    boundary = bottom
    value = 0.0
  [../]
  [./back_z]
    type = DirichletBC
    variable = disp_z
    boundary = 'back'
    value = 0.0
  [../]
[]

[Materials]
  [./Copper]
    type = GBEvolution
    block = 0
    T = 500 # K
    wGB = 15 # nm
    GBmob0 = 2.5e-6 # m^4/(Js) from Schoenfelder 1997
    Q = 0.23 # Migration energy in eV
    GBenergy = 0.708 # GB energy in J/m^2
    # length_scale = 1.0e-6
    # time_scale = 1.0
  [../]
  [./ElasticityTensor]
    type = ComputePolycrystalElasticityTensor
    grain_tracker = grain_tracker
  [../]
  [./strain]
    type = ComputeSmallStrain
    block = 0
    displacements = 'disp_x disp_y disp_z'
  [../]
  [./stress]
    type = ComputeLinearElasticStress
    block = 0
  [../]
[]

[Postprocessors]
  [./ngrains]
    type = FeatureFloodCount
    variable = bnds
    threshold = 0.7
  [../]
  [./dofs]
    type = NumDOFs
  [../]
  [./dt]
    type = TimestepSize
  [../]
  [./run_time]
    type = PerfGraphData
    section_name = "Root"
    data_type = total
  [../]
  [./avg_grain_volumes]
    type = AverageGrainVolume
    feature_counter = grain_tracker
    execute_on = 'initial timestep_end'
  [../]
[]

[VectorPostprocessors]
  [./grain_volumes]
    type = FeatureVolumeVectorPostprocessor
    flood_counter = grain_tracker
    execute_on = 'initial timestep_end'
  #  output_centroids = true
  [../]
[]

[Preconditioning]
  [./SMP]
    type = SMP
    coupled_groups = 'disp_x,disp_y,disp_z'
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
  end_time = 5000.0

  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 1.5
    growth_factor = 1.2
    cutback_factor = 0.8
    optimal_iterations = 8
  [../]

  [./Adaptivity]
    initial_adaptivity = 2
    refine_fraction = 0.8
    coarsen_fraction = 0.05
    max_h_level = 3
  [../]
[]

[Outputs]
  file_base = ./${my_filename}/out_${my_filename}
 # csv = true
  [./csv]
    type = CSV
    interval = 4
  [../]
  [./exodus]
    type = Exodus
    interval = 2
  [../]
  [./console]
    type = Console
    max_rows = 20 # Will print the 20 most recent postprocessor values to the screen
  [../]
  [out]
    type = Checkpoint
    interval = 8
    num_files = 3
  []
  [pgraph]
    type = PerfGraphOutput
    execute_on = 'initial final'  # Default is "final"
    level = 2                     # Default is 1
    heaviest_branch = true        # Default is false
    heaviest_sections = 7         # Default is 0
  []
[]


# =0.5*(elastic_strain11*stress11+2*elastic_strain12*stress12+elastic_strain22*stress22)