#
# KKS simple example in the split form
# step03: Elastic energy coupled with Allen-Cahn

cs_eq='1.0'
cl_eq ='0.0'
interval_exodus = '25'
num_steps='200'
dt='0.1'
width='4.0'
radius='15.0'
length='64'

misfit=0.03
c11=40.6
c12=30
c44=27


[GlobalParams]
    int_width = ${width}
    radius    = ${radius}
    x1        = ${fparse ${length} / 2.0}
    y1        = ${fparse ${length} / 2.0}
    end_point ='0   ${fparse ${length} / 2.0 } 0 '
    start_point   ='${fparse ${length} / 2.0 } ${fparse ${length} / 2.0 } 0' 
    displacements = 'disp_x disp_y'
    enable_jit = true
    use_legacy_material_output = false
[]

[Mesh]
        [./gen]
        type = GeneratedMeshGenerator
  dim = 2
  nx = ${length}
  ny = ${length}
  xmax = ${length}
  ymax = ${length}
        [../]

  [./cnode]
    type = ExtraNodesetGenerator
    input = gen
    coord = '${fparse ${length} / 2.0} ${fparse ${length} / 2.0}'
    new_boundary = 100
    tolerance = 0.1
  [../]
[]

[AuxVariables]
  [./Fglobal]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_vm]
    order = CONSTANT
    family= MONOMIAL
  [../]
[]

[Variables]
  # order parameter
  [./eta]
    order = FIRST
    family = LAGRANGE
  [../]

  # solute concentration
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
    initial_condition = '${cl_eq}'
  [../]
  # Solid phase solute concentration
  [./cs]
    order = FIRST
    family = LAGRANGE
    initial_condition = '${cs_eq}'
  [../]
[]


[Materials]
  # Free energy of the liquid
  [./fl]
    type = DerivativeParsedMaterial
    f_name = fl
    args = 'cl'
    function = '(${cl_eq}-cl)^2'
  [../]

  # Free energy of the solid
  [./fs]
    type = DerivativeParsedMaterial
    f_name = fs
    args = 'cs'
    function = '(${cs_eq}-cs)^2'
  [../]

  # h(eta)
  [./h_eta]
    type = SwitchingFunctionMaterial
    h_order = HIGH
    eta = eta
  [../]

  # g(eta)
  [./g_eta]
    type = BarrierFunctionMaterial
    g_order = SIMPLE
    eta = eta
  [../]

  # constant properties
  [./constants]
    type = GenericConstantMaterial
    prop_names  = 'M   L   eps_sq'
    prop_values = '1.0 1.0 1.0  '
  [../]

#Mechanics parameeters  
  [./Stiffness_matrix]
    type = ComputeElasticityTensor
    C_ijkl = '${c11} ${c12} ${c12} ${c11} ${c12} ${c11} ${c44} ${c44} ${c44}'
    fill_method = symmetric9
  [../]

  [./stress]
    type = ComputeLinearElasticStress
  [../]

  [./eigen_strain]
    type = ComputeVariableEigenstrain
    eigen_base = '${misfit} ${misfit} 0 0 0 0'
    prefactor = h
    args = eta
    eigenstrain_name = 'eigenstrain_ppt'
  [../]

  [./elastic_free_energy_p]
    type = ElasticEnergyMaterial
    f_name = f_el_mat
    args = 'eta'
    derivative_order = 3
  [../]
[]

[Kernels]
  active = 'PhaseConc ChemPotSolute CHBulk ACBulkF ACBulkC ACBulk_el ACInterface dcdt detadt ckernel'

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
  [./ACBulk_el] #This adds df_el/deta for strain interpolation
    type = AllenCahn
    variable = eta
    f_name = f_el_mat
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
  [./compute_stress_vm]
    variable = stress_vm
    type = ParsedAux
    function = 'sqrt(stress_xx^2+stress_yy^2-stress_xx*stress_yy+3*stress_xy^2)/sqrt(2)'
    args = 'stress_xx stress_yy stress_xy'
    execute_on = 'INITIAL TIMESTEP_END'
  [../]
[]

[Executioner]
  type = Transient
  solve_type = 'PJFNK'

  petsc_options_iname = '-pc_type -sub_pc_type -sub_pc_factor_shift_type'
  petsc_options_value = 'asm      ilu          nonzero'

  num_steps = ${num_steps}
  dt = ${dt}

  l_max_its = 100
  nl_max_its = 100
  l_tol = 1e-3
  l_abs_tol=1e-4
  nl_rel_tol=1e-8

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



[Outputs]
        [./exo]
                type = Exodus
                interval=${interval_exodus}
                elemental_as_nodal = 'true'
        [../]
[]


[ICs]
        [./eta]
                type = SmoothCircleIC
                invalue = 1.0
                outvalue = 0.0
                variable = eta
        [../]
        [./conc]
                type = SmoothCircleIC
                invalue = ${cs_eq}
                outvalue = ${cl_eq}
                variable = c
        [../]
[]


[BCs]
        [./Periodic]
                [./all]
                        variable       = 'c eta w'
                        auto_direction = 'x y'
                [../]
                [./mechanics]
                        variable = 'disp_x disp_y'
                        auto_direction='x y'
                [../]
        [../]

# Stress free boundary condition        
  [./centerfix_x]
    type = DirichletBC
    boundary = 100
    variable = disp_x
    value = 0
  [../]
  # fix side point x coordinate to inhibit rotation
  [./angularfix_y]
    type = DirichletBC
    boundary = bottom
    variable = disp_y
    value = 0
  [../]
[]



[Modules]
        [./TensorMechanics]
                [./Master]
                        [./all]
  strain = SMALL
  add_variables = true
    displacements = 'disp_x disp_y'
  eigenstrain_names = 'eigenstrain_ppt'
  generate_output = 'stress_xx stress_yy stress_xy elastic_strain_xx elastic_strain_yy elastic_strain_xy '
                        [../]
                [../]
        [../]
[]
