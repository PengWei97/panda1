my_filename = 'GG_exmaples200_elastic_3D' 
my_function = 'if(t<4,0.1*t,0.4+0.002*sin(10*pi*t))' # 0.001s^{-1}

[Mesh]
  type = GeneratedMesh
  dim = 3
  nx = 40
  ny = 40
  nz = 1
  xmin = 0
  xmax = 2000
  ymin = 0
  ymax = 2000
  zmin = 0
  zmax = 25
  elem_type = HEX8
  uniform_refine = 2
[]

[GlobalParams]
  op_num = 8
  var_name_base = gr
  grain_num = 200
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
    columnar_3D = true
    coloring_algorithm = jp
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
[]

[BCs]
  [./Periodic]
    [./All]
      auto_direction = 'x'
      variable = 'gr0 gr1 gr2 gr3 gr4 gr5 gr6 gr7'
    [../]
  [../]
  [./top_displacement]
    type = FunctionDirichletBC
    variable = disp_y
    boundary = top
    function = ${my_function}
  [../]
  [./x_anchor]
    type = DirichletBC
    variable = disp_x
    boundary = left
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
    # type = GBEvolution
    # T = 500 # K
    # wGB = 60 # nm
    # GBmob0 = 2.5e-6 #m^4/(Js) from Schoenfelder 1997
    # Q = 0.23 #Migration energy in eV
    # GBenergy = 0.708 #GB energy in J/m^2

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
  # active = ''
  [./ngrains]
    type = FeatureFloodCount
    variable = bnds
    threshold = 0.7
  [../]
  [./run_time]
    type = PerfGraphData
    section_name = "Root"
    data_type = total
  [../]
  [./dt]
    type = TimestepSize
  [../]
  [./avg_grain_volumes]
    type = AverageGrainVolume
    feature_counter = grain_tracker
    execute_on = 'initial timestep_end'
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
  scheme = 'bdf2'
  solve_type = 'PJFNK'

  petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart'
  petsc_options_value = 'hypre boomeramg 31'

  l_tol = 1.0e-5
  l_max_its = 15
  nl_max_its = 20
  nl_rel_tol = 1.0e-10
  start_time = 0.0
  dt = 5000.0

  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.01
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
  # exodus = true
  [./csv]
    type = CSV
    # interval = 1
  [../]
  [./exodus]
    type = Nemesis
    # interval = 2
  [../]
  [pgraph]
    type = PerfGraphOutput
    execute_on = 'initial final'  # Default is "final"
    level = 2                     # Default is 1
    heaviest_branch = true        # Default is false
    heaviest_sections = 7         # Default is 0
  []
[]
