./tutorials/darcy_thermo_mech/step09_mechanics/problems/step9.i:    function = 'if(t<0,350+50*t,350)'
./tutorials/darcy_thermo_mech/step09_mechanics/problems/step9.i:    function = 'if(t<0,0.1,0.25)'
./tutorials/darcy_thermo_mech/step07_adaptivity/problems/step7a_coarse.i:    function = 'if(t<0,350+50*t,350)'
./tutorials/darcy_thermo_mech/step07_adaptivity/problems/step7a_coarse.i:    function = 'if(t<0,0.1,0.25)'
./tutorials/darcy_thermo_mech/step07_adaptivity/problems/step7b_fine.i:    function = 'if(t<0,350+50*t,350)'
./tutorials/darcy_thermo_mech/step07_adaptivity/problems/step7b_fine.i:    function = 'if(t<0,0.1,0.25)'
./tutorials/darcy_thermo_mech/step07_adaptivity/problems/step7d_adapt_blocks.i:    function = 'if(t<0,350+50*t,350)'
./tutorials/darcy_thermo_mech/step07_adaptivity/problems/step7d_adapt_blocks.i:    function = 'if(t<0,0.1,0.25)'
./tutorials/darcy_thermo_mech/step07_adaptivity/problems/step7c_adapt.i:    function = 'if(t<0,350+50*t,350)'
./tutorials/darcy_thermo_mech/step07_adaptivity/problems/step7c_adapt.i:    function = 'if(t<0,0.1,0.25)'
./tutorials/darcy_thermo_mech/step08_postprocessors/problems/step8.i:    function = 'if(t<0,350+50*t,350)'
./tutorials/darcy_thermo_mech/step08_postprocessors/problems/step8.i:    function = 'if(t<0,0.1,0.25)'
./tutorials/darcy_thermo_mech/step06_coupled_darcy_heat_conduction/problems/step6a_coupled.i:    function = 'if(t<0,350+50*t,350)'
./tutorials/darcy_thermo_mech/step06_coupled_darcy_heat_conduction/problems/step6a_coupled.i:    function = 'if(t<0,0.1,0.25)'
./tutorials/darcy_thermo_mech/step06_coupled_darcy_heat_conduction/problems/step6b_transient_inflow.i:    function = 'if(t<0,350+50*t,350)'
./tutorials/darcy_thermo_mech/step06_coupled_darcy_heat_conduction/problems/step6b_transient_inflow.i:    function = 'if(t<0,0.1,(2*pi/(0.466*pi))/16)' # dt to always hit the peaks of sine/cosine BC
./tutorials/darcy_thermo_mech/step11_action/problems/step11.i:    function = 'if(t<0,350+50*t,350)'
./tutorials/darcy_thermo_mech/step11_action/problems/step11.i:    function = 'if(t<0,0.1,0.25)'
./tutorials/darcy_thermo_mech/step10_multiapps/problems/step10.i:    function = 'if(t<0,350+50*t,350)'
./tutorials/darcy_thermo_mech/step10_multiapps/problems/step10.i:    function = 'if(t<0,0.1,0.25)'
./test/tests/functions/piecewise_multilinear/twoD_const.i:            ix := if(x < 0.5, 0, if(x < 1, 1, 2));
./test/tests/functions/piecewise_multilinear/twoD_const.i:            iy := if(y > 0, 2, if(y > -0.5, 1, 0));
./test/tests/kernels/conservative_advection/none_in_none_out.i:    function = 'if(x<5,x,10-x)'
./test/tests/kernels/conservative_advection/no_upwinding_2D.i:    function = 'if(x<0.2,if(y<0.2,1,0),0)'
./test/tests/kernels/conservative_advection/full_upwinding_2D.i:    function = 'if(x<0.2,if(y<0.2,1,0),0)'
./test/tests/kernels/conservative_advection/none_in_all_out.i:    function = 'if(x<5,x,10-x)'
./test/tests/kernels/vector_fe/coupled_vector_gradient.i:    value = 'if(t < 1, 0, t - 1)'
./test/tests/kernels/vector_fe/coupled_vector_gradient.i:    value = 'if(t < 1, 0, t - 1)'
./test/tests/kernels/vector_fe/coupled_vector_gradient.i:    value = 'if(t < 2, 0, t - 2)'
./test/tests/kernels/vector_fe/coupled_vector_gradient.i:    value = 'if(t < 2, 0, t - 2)'
./test/tests/auxkernels/grad_component/grad_component_monomial.i:    function = 'if(x>0.5,if(x<1.5,2*x,3),0)'
./test/tests/userobjects/solution_user_object/discontinuous_value_solution_uo_p1.i:    value = 'if(x<0.5,3,5)'
./test/tests/userobjects/solution_user_object/discontinuous_value_solution_uo_p1.i:    value = 'if(x<0.5,x,2*x-0.5)'
./test/tests/userobjects/postprocessor_spatial_user_object/sub.i:    value = 'if(a < 0.8625, 1, 0)'
./test/tests/multiapps/picard_multilevel/multilevel_dt_rejection/picard_sub.i:    value = 'if(t < 2.5, 1, 1 / t)'
./test/tests/multiapps/grid-sequencing/vi-fine-alone.i:    function = 'if(x<5,-1,1)'
./test/tests/multiapps/grid-sequencing/vi-coarse.i:    function = 'if(x<5,-1,1)'
./test/tests/multiapps/grid-sequencing/vi-fine.i:    function = 'if(x<5,-1,1)'
./test/tests/multiapps/grid-sequencing/vi-coarser.i:    function = 'if(x<5,-1,1)'
./test/tests/dirackernels/function_dirac_source/function_dirac_source.i:    value = 'if(t < 1.0001, 1, 0)'
./test/tests/postprocessors/element_l1_error/element_l1_error.i:    value = 'if(x<5,5,6)'
./test/tests/postprocessors/element_l1_error/element_l1_error.i:    value = 'if(x<5,3,2)'
./test/tests/postprocessors/time_extreme_value/time_extreme_value.i:    function = 'if(t<1.0,t,1.0)'
./test/tests/postprocessors/time_extreme_value/time_extreme_value.i:    function = 'if(t<1.0,2.0-t,1.0)'
./test/tests/postprocessors/find_value_on_line/findvalueonline.i:      function = if(x<1,1-x,0)
./test/tests/utils/mathutils/smootherstep.i:             if(x < 0.2, 0, if(x > 0.8, 1, val))'
./test/tests/utils/mathutils/smootherstep.i:             if(x < 0.2, 0, if(x > 0.8, 0, val / (0.8 - 0.2)))'
./test/tests/fviks/one-var-diffusion/no-ik.i:    value = 'if(x<1, 1 - x/3, 4/3 - 2*x/3)'
./test/tests/fviks/one-var-diffusion/test.i:    value = 'if(x<1, 1 - x/3, 4/3 - 2*x/3)'
./test/tests/controls/bool_function_control/bool_function_control.i:    value = 'if(t<0.3, 1, 0)'
./test/tests/nodalkernels/constraint_enforcement/ad-upper-and-lower-bound.i:    function = 'if(x<5,-1,1)'
./test/tests/nodalkernels/constraint_enforcement/upper-and-lower-bound.i:    function = 'if(x<5,-1,1)'
./test/tests/nodalkernels/constraint_enforcement/vi-bounding.i:    function = 'if(x<5,-1,1)'
./modules/navier_stokes/test/tests/scalar_adr/supg/tauOpt.i:    value = 'if(x < 6, 1 - .25 * x, if(x < 8, -2 + .25 * x, 0))'
./modules/navier_stokes/test/tests/finite_volume/cns/straight_channel_porosity_step/dc.i:    value = 'if(y < 2.8, 1,
./modules/navier_stokes/test/tests/finite_volume/cns/straight_channel_porosity_step/dc.i:             if(y < 3.2, 1 - .5 / .4 * (y - 2.8),
./modules/navier_stokes/test/tests/finite_volume/cns/straight_channel_porosity_step/dc.i:             if(y < 6.8, .5,
./modules/navier_stokes/test/tests/finite_volume/cns/straight_channel_porosity_step/dc.i:             if(y < 7.2, .5 - .25 / .4 * (y - 6.8),
./modules/navier_stokes/test/tests/finite_volume/cns/straight_channel_porosity_step/dc.i:             if(y < 10.8, .25,
./modules/navier_stokes/test/tests/finite_volume/cns/straight_channel_porosity_step/dc.i:             if(y < 11.2, .25 + .25 / .4 * (y - 10.8),
./modules/navier_stokes/test/tests/finite_volume/cns/straight_channel_porosity_step/dc.i:             if(y < 14.8, .5,
./modules/navier_stokes/test/tests/finite_volume/cns/straight_channel_porosity_step/dc.i:             if(y < 15.2, .5 + .5 / .4 * (y - 14.8),
./modules/navier_stokes/test/tests/finite_volume/cns/straight_channel_porosity_step/implicit-euler-basic-kt-primitive.i:    value = 'if(x < 2, 1,
./modules/navier_stokes/test/tests/finite_volume/cns/straight_channel_porosity_step/implicit-euler-basic-kt-primitive.i:             if(x < 4, 1 - .5 / 2 * (x - 2),
./modules/navier_stokes/test/tests/finite_volume/cns/straight_channel_porosity_step/implicit-euler-basic-kt-primitive.i:             if(x < 6, .5,
./modules/navier_stokes/test/tests/finite_volume/cns/straight_channel_porosity_step/implicit-euler-basic-kt-primitive.i:             if(x < 8, .5 - .25 / 2 * (x - 6),
./modules/navier_stokes/test/tests/finite_volume/cns/straight_channel_porosity_step/implicit-euler-basic-kt-primitive.i:             if(x < 10, .25,
./modules/navier_stokes/test/tests/finite_volume/cns/straight_channel_porosity_step/implicit-euler-basic-kt-primitive.i:             if(x < 12, .25 + .25 / 2 * (x - 10),
./modules/navier_stokes/test/tests/finite_volume/cns/straight_channel_porosity_step/implicit-euler-basic-kt-primitive.i:             if(x < 14, .5,
./modules/navier_stokes/test/tests/finite_volume/cns/straight_channel_porosity_step/implicit-euler-basic-kt-primitive.i:             if(x < 16, .5 + .5 / 2 * (x - 14),
./modules/navier_stokes/test/tests/finite_volume/cns/straight_channel_porosity_step/rotated-2d-bkt-function-porosity.i:    value = 'if(y < 2.8, 1,
./modules/navier_stokes/test/tests/finite_volume/cns/straight_channel_porosity_step/rotated-2d-bkt-function-porosity.i:             if(y < 3.2, 1 - .5 / .4 * (y - 2.8),
./modules/navier_stokes/test/tests/finite_volume/cns/straight_channel_porosity_step/rotated-2d-bkt-function-porosity.i:             if(y < 6.8, .5,
./modules/navier_stokes/test/tests/finite_volume/cns/straight_channel_porosity_step/rotated-2d-bkt-function-porosity.i:             if(y < 7.2, .5 - .25 / .4 * (y - 6.8),
./modules/navier_stokes/test/tests/finite_volume/cns/straight_channel_porosity_step/rotated-2d-bkt-function-porosity.i:             if(y < 10.8, .25,
./modules/navier_stokes/test/tests/finite_volume/cns/straight_channel_porosity_step/rotated-2d-bkt-function-porosity.i:             if(y < 11.2, .25 + .25 / .4 * (y - 10.8),
./modules/navier_stokes/test/tests/finite_volume/cns/straight_channel_porosity_step/rotated-2d-bkt-function-porosity.i:             if(y < 14.8, .5,
./modules/navier_stokes/test/tests/finite_volume/cns/straight_channel_porosity_step/rotated-2d-bkt-function-porosity.i:             if(y < 15.2, .5 + .5 / .4 * (y - 14.8),
./modules/navier_stokes/test/tests/finite_volume/cns/straight_channel_porosity_step/rotated-2d-bkt-function-porosity-mixed.i:    value = 'if(y < 2.8, 1,
./modules/navier_stokes/test/tests/finite_volume/cns/straight_channel_porosity_step/rotated-2d-bkt-function-porosity-mixed.i:             if(y < 3.2, 1 - .5 / .4 * (y - 2.8),
./modules/navier_stokes/test/tests/finite_volume/cns/straight_channel_porosity_step/rotated-2d-bkt-function-porosity-mixed.i:             if(y < 6.8, .5,
./modules/navier_stokes/test/tests/finite_volume/cns/straight_channel_porosity_step/rotated-2d-bkt-function-porosity-mixed.i:             if(y < 7.2, .5 - .25 / .4 * (y - 6.8),
./modules/navier_stokes/test/tests/finite_volume/cns/straight_channel_porosity_step/rotated-2d-bkt-function-porosity-mixed.i:             if(y < 10.8, .25,
./modules/navier_stokes/test/tests/finite_volume/cns/straight_channel_porosity_step/rotated-2d-bkt-function-porosity-mixed.i:             if(y < 11.2, .25 + .25 / .4 * (y - 10.8),
./modules/navier_stokes/test/tests/finite_volume/cns/straight_channel_porosity_step/rotated-2d-bkt-function-porosity-mixed.i:             if(y < 14.8, .5,
./modules/navier_stokes/test/tests/finite_volume/cns/straight_channel_porosity_step/rotated-2d-bkt-function-porosity-mixed.i:             if(y < 15.2, .5 + .5 / .4 * (y - 14.8),
./modules/navier_stokes/test/tests/finite_volume/cns/straight_channel_porosity_step/hllc.i:    value = 'if(x < 2, 1,
./modules/navier_stokes/test/tests/finite_volume/cns/straight_channel_porosity_step/hllc.i:             if(x < 4, 1 - .5 / 2 * (x - 2),
./modules/navier_stokes/test/tests/finite_volume/cns/straight_channel_porosity_step/hllc.i:             if(x < 6, .5,
./modules/navier_stokes/test/tests/finite_volume/cns/straight_channel_porosity_step/hllc.i:             if(x < 8, .5 - .25 / 2 * (x - 6),
./modules/navier_stokes/test/tests/finite_volume/cns/straight_channel_porosity_step/hllc.i:             if(x < 10, .25,
./modules/navier_stokes/test/tests/finite_volume/cns/straight_channel_porosity_step/hllc.i:             if(x < 12, .25 + .25 / 2 * (x - 10),
./modules/navier_stokes/test/tests/finite_volume/cns/straight_channel_porosity_step/hllc.i:             if(x < 14, .5,
./modules/navier_stokes/test/tests/finite_volume/cns/straight_channel_porosity_step/hllc.i:             if(x < 16, .5 + .5 / 2 * (x - 14),
./modules/phase_field/tutorials/spinodal_decomposition/s5_energycurve.i:    function = if(c>0.6,0.0016,0)
./modules/phase_field/tutorials/spinodal_decomposition/s4_mobility.i:    function = if(c>0.6,0.0016,0)
./modules/phase_field/test/tests/MultiPhase/asymmetriccrosstermbarrierfunction.i:    value = 'il:=x-7; ir:=2-x; if(x<1, 1,
./modules/phase_field/test/tests/MultiPhase/asymmetriccrosstermbarrierfunction.i:                               if(x<2, 0.5-0.5*cos(ir*pi),
./modules/phase_field/test/tests/MultiPhase/asymmetriccrosstermbarrierfunction.i:                               if(x<7, 0,
./modules/phase_field/test/tests/MultiPhase/asymmetriccrosstermbarrierfunction.i:                               if(x<8, 0.5-0.5*cos(il*pi),
./modules/phase_field/test/tests/MultiPhase/asymmetriccrosstermbarrierfunction.i:    value = 'il:=x-1; ir:=5-x; if(x<1, 0,
./modules/phase_field/test/tests/MultiPhase/asymmetriccrosstermbarrierfunction.i:                               if(x<2, 0.5-0.5*cos(il*pi),
./modules/phase_field/test/tests/MultiPhase/asymmetriccrosstermbarrierfunction.i:                               if(x<4, 1,
./modules/phase_field/test/tests/MultiPhase/asymmetriccrosstermbarrierfunction.i:                               if(x<5, 0.5-0.5*cos(ir*pi),
./modules/phase_field/test/tests/MultiPhase/asymmetriccrosstermbarrierfunction.i:    value = 'il:=x-4; ir:=8-x; if(x<4, 0,
./modules/phase_field/test/tests/MultiPhase/asymmetriccrosstermbarrierfunction.i:                               if(x<5, 0.5-0.5*cos(il*pi),
./modules/phase_field/test/tests/MultiPhase/asymmetriccrosstermbarrierfunction.i:                               if(x<7, 1,
./modules/phase_field/test/tests/MultiPhase/asymmetriccrosstermbarrierfunction.i:                               if(x<8, 0.5-0.5*cos(ir*pi),
./modules/phase_field/test/tests/MultiPhase/crosstermbarrierfunction.i:    value = 'il:=x-7; ir:=2-x; if(x<1, 1,
./modules/phase_field/test/tests/MultiPhase/crosstermbarrierfunction.i:                               if(x<2, 0.5-0.5*cos(ir*pi),
./modules/phase_field/test/tests/MultiPhase/crosstermbarrierfunction.i:                               if(x<7, 0,
./modules/phase_field/test/tests/MultiPhase/crosstermbarrierfunction.i:                               if(x<8, 0.5-0.5*cos(il*pi),
./modules/phase_field/test/tests/MultiPhase/crosstermbarrierfunction.i:    value = 'il:=x-1; ir:=5-x; if(x<1, 0,
./modules/phase_field/test/tests/MultiPhase/crosstermbarrierfunction.i:                               if(x<2, 0.5-0.5*cos(il*pi),
./modules/phase_field/test/tests/MultiPhase/crosstermbarrierfunction.i:                               if(x<4, 1,
./modules/phase_field/test/tests/MultiPhase/crosstermbarrierfunction.i:                               if(x<5, 0.5-0.5*cos(ir*pi),
./modules/phase_field/test/tests/MultiPhase/crosstermbarrierfunction.i:    value = 'il:=x-4; ir:=8-x; if(x<4, 0,
./modules/phase_field/test/tests/MultiPhase/crosstermbarrierfunction.i:                               if(x<5, 0.5-0.5*cos(il*pi),
./modules/phase_field/test/tests/MultiPhase/crosstermbarrierfunction.i:                               if(x<7, 1,
./modules/phase_field/test/tests/MultiPhase/crosstermbarrierfunction.i:                               if(x<8, 0.5-0.5*cos(ir*pi),
./modules/phase_field/test/tests/MaskedBodyForce/MaskedBodyForce_test.i:    function = if(c>0.5,0,1)
./modules/phase_field/test/tests/SoretDiffusion/split_temp.i:    function = 'if(c>0.7,1e-8,4e-8)'
./modules/phase_field/test/tests/SoretDiffusion/direct_temp.i:    function = 'if(c>0.7,1e-8,4e-8)'
./modules/phase_field/test/tests/misc/equal_gradient_lagrange.i:      function = 'r:=sqrt((x-0.4)^2+(y-0.5)^2);if(r<0.05,5,1)'
./modules/phase_field/test/tests/misc/interface_flux.i:      function = 'r:=sqrt((x-0.4)^2+(y-0.5)^2);if(r<0.05,5,1)'
./modules/phase_field/test/tests/misc/interface_flux.i:      function = 'r:=sqrt((x-0.7)^2+(y-0.5)^2);if(r<0.05,5,1)'
./modules/phase_field/test/tests/misc/interface_grad.i:      function = 'r:=sqrt((x-0.4)^2+(y-0.5)^2);if(r<0.05,5,1)'
./modules/phase_field/test/tests/mobility_derivative/AC_mobility_derivative_test.i:    function = 'if(op<0, 0.01, if(op>1, 0.01, 1*op^2*(1-op)^2+0.01))'
./modules/phase_field/test/tests/mobility_derivative/AC_mobility_derivative_coupled_test.i:    function = 'l:=0.1+1*(v+op)^2; if(l<0.01, 0.01, l)'
./modules/phase_field/test/tests/mobility_derivative/mobility_derivative_direct_coupled_test.i:    function = if(d>0.001,d,0.001)*if(c<0,0.5,if(c>1,0.5,1-0.5*c^2))
./modules/phase_field/test/tests/mobility_derivative/mobility_derivative_split_coupled_test.i:    function = 'if(d>0.001,d,0.001)*(1-0.5*c^2)'
./modules/phase_field/test/tests/mobility_derivative/mobility_derivative_direct_test.i:    function = 'if(c<-1,0.1,if(c>1,0.1,1-.9*c^2))'
./modules/phase_field/test/tests/grain_boundary_area/diagonal.i:      function = 'd:=(x-y)*80;if(d<pi&d>-pi,sin(d/2)/2+0.5,if(d<0,0,1))'
./modules/phase_field/test/tests/grain_boundary_area/diagonal.i:      function = 'd:=(x-y)*80;1-if(d<pi&d>-pi,sin(d/2)/2+0.5,if(d<0,0,1))'
./modules/phase_field/examples/multiphase/DerivativeMultiPhaseMaterial.i:      function = 'r:=sqrt(x^2+y^2);if(r<=4,1,0)'
./modules/phase_field/examples/multiphase/DerivativeMultiPhaseMaterial.i:      function = 'r:=sqrt(x^2+y^2);if(r>4&r<=7,1,0)'
./modules/phase_field/examples/multiphase/DerivativeMultiPhaseMaterial.i:      function = 'r:=sqrt(x^2+y^2);if(r>7,1,0)'
./modules/phase_field/examples/nucleation/refine.i:    function = 'if(c<0.21,c*1e-8,0)'
./modules/phase_field/examples/interfacekernels/interface_gradient.i:      function = 'r:=sqrt((x-0.4)^2+(y-0.5)^2);if(r<0.05,5,1)'
./modules/phase_field/examples/interfacekernels/interface_fluxbc.i:      function = 'r:=sqrt((x-0.4)^2+(y-0.5)^2);if(r<0.05,5,1)'
./modules/phase_field/examples/interfacekernels/interface_fluxbc.i:      function = 'r:=sqrt((x-0.7)^2+(y-0.5)^2);if(r<0.05,5,1)'
./modules/phase_field/examples/interfacekernels/interface_fluxbc.i:      function = 'r:=sqrt((x-0.4)^2+(y-0.5)^2);if(r<0.05,5,1)'
./modules/phase_field/examples/interfacekernels/interface_fluxbc.i:      function = 'r:=sqrt((x-0.7)^2+(y-0.5)^2);if(r<0.05,5,1)'
./modules/contact/test/tests/non-singular-frictional-mortar/frictional-mortar.i:    value = 'if(t<1.0,${vx}*t-${offset},${vx}-${offset})'
./modules/contact/test/tests/non-singular-frictional-mortar/frictional-mortar.i:    value = 'if(t<1.0,${offset},${vy}*(t-1.0)+${offset})'
./modules/contact/test/tests/fieldsplit/frictional_mortar_FS.i:    value = 'if(t<0.5,${vx}*t-${offset},${vx}-${offset})'
./modules/contact/test/tests/fieldsplit/frictional_mortar_FS.i:    value = 'if(t<0.5,${offset},${vy}*(t-0.5)+${offset})'
./modules/contact/examples/3d_berkovich/indenter_berkovich_friction.i:    value = 'if(t < 1.5, -t, t-3.0)'
./modules/rdg/test/tests/advection_1d/rdgP0.i:    function = 'if(x<0.1,0,if(x>0.3,0,1))'
./modules/heat_conduction/test/tests/sideset_heat_transfer/gap_thermal_ktemp_1D.i:            A := if(x < 1, -0.5, -0.25);
./modules/heat_conduction/test/tests/sideset_heat_transfer/gap_thermal_ktemp_1D.i:            B := if(x < 1, -0.293209850655001, 0.0545267662299068);
./modules/heat_conduction/test/tests/sideset_heat_transfer/gap_thermal_ktemp_1D.i:            C := if(x < 1, 300.206790149345, 300.19547323377);
./modules/heat_conduction/test/tests/sideset_heat_transfer/gap_thermal_1D.i:            A := if(x < 1, -0.5, -0.25);
./modules/heat_conduction/test/tests/sideset_heat_transfer/gap_thermal_1D.i:            B := if(x < 1, -0.293209850655001, 0.0545267662299068);
./modules/heat_conduction/test/tests/sideset_heat_transfer/gap_thermal_1D.i:            C := if(x < 1, 300.206790149345, 300.19547323377);
./modules/combined/test/tests/poro_mechanics/terzaghi.i:    function = if(0.5*t<0.1,0.5*t,0.1)
./modules/combined/test/tests/poro_mechanics/mandel.i:    function = if(0.15*t<0.01,0.15*t,0.01)
./modules/combined/test/tests/phase_field_fracture/crack2d_vi_solver.i:    value = 'if(x<0.5 & y < 0.55 & y > 0.45,1, 0)'
./modules/combined/examples/effective_properties/effective_th_cond.i:    function = 'sk_b:= length_scale*k_b; sk_p2:= length_scale*k_p2; sk_int:= k_int*length_scale; if(phase2>0.1,if(phase2>0.95,sk_p2,sk_int),sk_b)'
./modules/combined/examples/phase_field-mechanics/interface_stress.i:    value = 'r:=sqrt(x^2+y^2+z^2); R:=(4.0-r)/2.0; if(R>1,1,if(R<0,0,3*R^2-2*R^3))'
./modules/combined/examples/phase_field-mechanics/Pattern1.i:    function = 'if(eta3>0.5,1,0)-if(eta2>0.5,1,0)'
./modules/combined/examples/mortar/mortar_gradient.i:    value = 'if(x>0.4 & x<0.6 & y>0.1 & y<0.3, 3+y, y)'
./modules/combined/examples/mortar/mortar_gradient.i:    value = 'if(x>0.4 & x<0.6 & y>0.1 & y<0.3, 3, 0)'
./modules/combined/examples/publications/rapid_dev/fig8.i:    function = 'if(eta3>0.5,1,0)-if(eta2>0.5,1,0)'
./modules/combined/examples/publications/rapid_dev/fig3.i:      #function = if(x>0,1,0)
./modules/combined/examples/geochem-porous_flow/geotes_2D/aquifer_un_quartz_geochemistry.i:    function = 'if(abs(x) = 70 & abs(y) = 40, 2.5, if(abs(x) = 70 | abs(y) = 40, 5, 10))'
./modules/combined/examples/geochem-porous_flow/forge/water_60_to_220degC.i:    function = 'if(t<1, 1, if(t<1.01, 0.01, 1))'
./modules/combined/examples/geochem-porous_flow/forge/water_60_to_220degC.i:    function = 'if(t>1, 16.8, 0)'
./modules/combined/examples/geochem-porous_flow/forge/water_60_to_220degC.i:    function = 'if(t>1, 1.8, 0)'
./modules/combined/examples/geochem-porous_flow/forge/water_60_to_220degC.i:    function = 'if(t>1, 10.4, 0)'
./modules/combined/examples/geochem-porous_flow/forge/water_60_to_220degC.i:    function = 'if(t>1, 30.0, 0)'
./modules/combined/examples/geochem-porous_flow/forge/water_60_to_220degC.i:    function = 'if(t>1, 0.48, 0)'
./modules/combined/examples/geochem-porous_flow/forge/water_60_to_220degC.i:    function = 'if(t>1, 0.52, 0)'
./modules/geochemistry/test/tests/kernels/dispersion_1.i:    function = 'if(x<=0.0, -1.0, 1.0)'
./modules/geochemistry/test/tests/kernels/dispersion_1.i:    value = 'xi := x / sqrt(4 * t * 0.3); expxi := exp(-xi * xi); if(x < 0.0, -1.0, if(x > 0.0, 1.0, 0.0)) * 2 / sqrt(pi) * sqrt(1 - expxi) * (sqrt(pi) / 2.0 + 31.0 * expxi / 200.0 - 341.0 * expxi * expxi / 8000.0)'
./modules/geochemistry/test/tests/kernels/time_deriv_2.i:   function = 'if(x > 0.5 & x < 1.5, x * x + 2.0 * 12.0 / (6.0 + x), x * x)'
./modules/geochemistry/test/tests/kernels/advection_1.i:    function = 'if(x<=0.25, 1, 0)'
./modules/geochemistry/test/tests/nodal_void_volume/nodal_void_volume.i:    function = 'if(x<4, 1, 2)'
./modules/geochemistry/test/tests/nodal_void_volume/nodal_void_volume_adaptive.i:    function = 'if(x<2,0,1)'
./modules/geochemistry/test/tests/nodal_void_volume/nodal_void_volume_adaptive.i:    function = 'if(x<4, 1, 2)'
./modules/geochemistry/test/tests/spatial_reactor/spatial_7.i:    function = 'if(t <= 1, 25, 25 + 18 * x)'
./modules/geochemistry/test/tests/time_dependent_reactions/seawater_evaporation_flow_through.i:    function = 'if(t<=1.0, 1.0, 2.0)' # initial "dump" then "flow_through"
./modules/geochemistry/test/tests/time_dependent_reactions/mixing.i:    function = 'if(t<=0, 1, 0)' # dump at start of first timestep
./modules/geochemistry/test/tests/time_dependent_reactions/mixing.i:    function = 'if(t<=0, 273, 4)' # during initialisation and dumping, T=273, while during adding T=temperature of reactants
./modules/geochemistry/test/tests/time_dependent_reactions/mixing.i:    function = 'if(t<=0, 0, 55.510000000000005)'
./modules/geochemistry/test/tests/time_dependent_reactions/mixing.i:    function = 'if(t<=0, 0, 3.643e-10)'
./modules/geochemistry/test/tests/time_dependent_reactions/mixing.i:    function = 'if(t<=0, 0, 8.831e-08)'
./modules/geochemistry/test/tests/time_dependent_reactions/mixing.i:    function = 'if(t<=0, 0, 0.0104)'
./modules/geochemistry/test/tests/time_dependent_reactions/mixing.i:    function = 'if(t<=0, 0, 0.559)'
./modules/geochemistry/test/tests/time_dependent_reactions/mixing.i:    function = 'if(t<=0, 0, 7.000000000000001e-09)'
./modules/geochemistry/test/tests/time_dependent_reactions/mixing.i:    function = 'if(t<=0, 0, 4.746e-15)'
./modules/geochemistry/test/tests/time_dependent_reactions/mixing.i:    function = 'if(t<=0, 0, 0.0002005)'
./modules/geochemistry/test/tests/time_dependent_reactions/mixing.i:    function = 'if(t<=0, 0, 0.002153)'
./modules/geochemistry/test/tests/time_dependent_reactions/mixing.i:    function = 'if(t<=0, 0, 0.010100000000000001)'
./modules/geochemistry/test/tests/time_dependent_reactions/mixing.i:    function = 'if(t<=0, 0, 0.054400000000000004)'
./modules/geochemistry/test/tests/time_dependent_reactions/mixing.i:    function = 'if(t<=0, 0, 6.79e-14)'
./modules/geochemistry/test/tests/time_dependent_reactions/mixing.i:    function = 'if(t<=0, 0, 0.48019999999999996)'
./modules/geochemistry/test/tests/time_dependent_reactions/mixing.i:    function = 'if(t<=0, 0, 0.000123)'
./modules/geochemistry/test/tests/time_dependent_reactions/mixing.i:    function = 'if(t<=0, 0, 0.0295)'
./modules/geochemistry/test/tests/time_dependent_reactions/mixing.i:    function = 'if(t<=0, 0, 0.00017)'
./modules/geochemistry/test/tests/time_dependent_reactions/mixing.i:    function = 'if(t<=0, 0, 3.8350000000000004e-05)'
./modules/geochemistry/test/tests/time_dependent_reactions/mixing.i:    function = 'if(t<=0, 0, 1e-08)'
./modules/geochemistry/test/tests/time_dependent_reactions/seawater_evaporation_no_flow_through.i:    function = 'if(t<=1.0, 1.0, 0.0)' # initial "dump" then "normal"
./modules/xfem/test/tests/moving_interface/phase_transition_2d.i:    function = 'if(x<5.01, 2, 1)'
./modules/xfem/test/tests/moving_interface/phase_transition_3d.i:    function = 'if(x<5.01, 2, 1)'
./modules/xfem/test/tests/moving_interface/ad_phase_transition_2d.i:    function = 'if(x<5.01, 2, 1)'
./modules/tensor_mechanics/test/tests/finite_strain_tensor_mechanics_tests/elastic_rotation.i:    value = 'if(t<=1.0, delta*t, (1.0+delta)*cos(pi/2*(t-t0)) - 1.0)'
./modules/tensor_mechanics/test/tests/finite_strain_tensor_mechanics_tests/elastic_rotation.i:    value = 'if(t<=1.0, 0.0, (1.0+delta)*sin(pi/2*(t-t0)))'
./modules/tensor_mechanics/test/tests/finite_strain_tensor_mechanics_tests/elastic_rotation.i:    value = 'if(t<=1.0, delta*t, (1.0+delta)*cos(pi/2.0*(t-t0)) - sin(pi/2.0*(t-t0)) - 1.0)'
./modules/tensor_mechanics/test/tests/finite_strain_tensor_mechanics_tests/elastic_rotation.i:    value = 'if(t<=1.0, 0.0, cos(pi/2.0*(t-t0)) + (1+delta)*sin(pi/2.0*(t-t0)) - 1.0)'
./modules/tensor_mechanics/test/tests/finite_strain_tensor_mechanics_tests/elastic_rotation.i:    value = 'if(t<=1.0, 0.0, -sin(pi/2.0*(t-t0)))'
./modules/tensor_mechanics/test/tests/finite_strain_tensor_mechanics_tests/elastic_rotation.i:    value = 'if(t<=1.0, 0.0, cos(pi/2.0*(t-t0)) - 1.0)'
./modules/tensor_mechanics/test/tests/capped_mohr_coulomb/random3.i:    value = 'if(a<1E-1,0,a)'
./modules/tensor_mechanics/test/tests/capped_mohr_coulomb/random3.i:    value = 'if(a<1E-1,0,a)'
./modules/tensor_mechanics/test/tests/capped_mohr_coulomb/random3.i:    value = 'if(a<1E-1,0,a)'
./modules/tensor_mechanics/test/tests/capped_mohr_coulomb/random3.i:    value = 'if(a<1E-1,0,a)'
./modules/tensor_mechanics/test/tests/capped_mohr_coulomb/random3.i:    value = 'if(a<1E-1,0,a)'
./modules/tensor_mechanics/test/tests/capped_mohr_coulomb/random3.i:    value = 'if(a<1E-1,0,a)'
./modules/tensor_mechanics/test/tests/capped_mohr_coulomb/random4.i:    value = 'if(a<1E-1,0,a)'
./modules/tensor_mechanics/test/tests/capped_mohr_coulomb/random4.i:    value = 'if(a<1E-1,0,a)'
./modules/tensor_mechanics/test/tests/capped_mohr_coulomb/random4.i:    value = 'if(a<1E-1,0,a)'
./modules/tensor_mechanics/test/tests/capped_mohr_coulomb/random4.i:    value = 'if(a<1E-1,0,a)'
./modules/tensor_mechanics/test/tests/capped_mohr_coulomb/random4.i:    value = 'if(a<1E-1,0,a)'
./modules/tensor_mechanics/test/tests/capped_mohr_coulomb/random4.i:    value = 'if(a<1E-1,0,a)'
./modules/tensor_mechanics/test/tests/capped_mohr_coulomb/random4.i:    value = 'if(a<1E-1,0,a)'
./modules/tensor_mechanics/test/tests/capped_mohr_coulomb/random4.i:    value = 'if(a<1E-1,0,a)'
./modules/tensor_mechanics/test/tests/capped_mohr_coulomb/random4.i:    value = 'if(a<1E-1,0,a)'
./modules/tensor_mechanics/test/tests/capped_mohr_coulomb/random4.i:    value = 'if(a<1E-1,0,a)'
./modules/tensor_mechanics/test/tests/capped_mohr_coulomb/random4.i:    value = 'if(a<1E-1,0,a)'
./modules/tensor_mechanics/test/tests/capped_mohr_coulomb/random4.i:    value = 'if(a<1E-1,0,a)'
./modules/tensor_mechanics/test/tests/capped_mohr_coulomb/random2.i:    value = 'if(a<1E-1,0,a)'
./modules/tensor_mechanics/test/tests/capped_mohr_coulomb/random2.i:    value = 'if(a<1E-1,0,a)'
./modules/tensor_mechanics/test/tests/capped_mohr_coulomb/random2.i:    value = 'if(a<1E-1,0,a)'
./modules/tensor_mechanics/test/tests/capped_mohr_coulomb/random1.i:    value = 'if(a<1E-1,0,a)'
./modules/tensor_mechanics/test/tests/capped_mohr_coulomb/random1.i:    value = 'if(a<1E-1,0,a)'
./modules/tensor_mechanics/test/tests/capped_mohr_coulomb/random1.i:    value = 'if(a<1E-1,0,a)'
./modules/tensor_mechanics/test/tests/tensile/random_smoothed.i:    value = 'if(a<1E-1,0,a)'
./modules/tensor_mechanics/test/tests/tensile/random_planar.i:    value = 'if(a<1E-1,0,a)'
./modules/tensor_mechanics/test/tests/tensile/random_planar.i:    value = 'if(a<1E-1,0,a)'
./modules/tensor_mechanics/test/tests/tensile/random_planar.i:    value = 'if(a<1E-1,0,a)'
./modules/tensor_mechanics/test/tests/tensile/random_update.i:    value = 'if(a<1E-1,0,a)'
./modules/tensor_mechanics/test/tests/tensile/random_update.i:    value = 'if(a<1E-1,0,a)'
./modules/tensor_mechanics/test/tests/tensile/random_update.i:    value = 'if(a<1E-1,0,a)'
./modules/tensor_mechanics/test/tests/multi/rock1.i:    value = 'if(a<1E-1,0,a)'
./modules/tensor_mechanics/test/tests/multi/rock1.i:    value = 'if(a<1E-1,0,a)'
./modules/tensor_mechanics/test/tests/multi/rock1.i:    value = 'if(a<1E-1,0,a)'
./modules/tensor_mechanics/test/tests/multi/rock1.i:    value = 'if(a<1E-1,0,a)'
./modules/tensor_mechanics/test/tests/smeared_cracking/cracking_rotation.i:    function = 'if(t<10,0,if(t>=100,1,1-cos((t-10)*pi/180)))'
./modules/tensor_mechanics/test/tests/smeared_cracking/cracking_rotation.i:    function = 'if(t<10,0,if(t>=100,1,sin((t-10)*pi/180)))'
./modules/tensor_mechanics/test/tests/smeared_cracking/cracking_rotation.i:    function = 'if(t<5,t*0.01,0.05-(t-5)*0.01)'
./modules/tensor_mechanics/test/tests/smeared_cracking/cracking_function.i:    value = 'if(x > 0.667, 1.1e6, 1.2e6)'
./modules/tensor_mechanics/test/tests/smeared_cracking/cracking_xyz.i:    value = 'if(t < 2, 0.00175, 0)'
./modules/tensor_mechanics/test/tests/drucker_prager/random_hyperbolic.i:    value = 'if(a<1E-3,0,a)'
./modules/tensor_mechanics/test/tests/rom_stress_update/ADlower_limit.i:             clamped := if(val <= -1, -0.99999, if(val >= 1, 0.99999, val));
./modules/tensor_mechanics/test/tests/rom_stress_update/ADlower_limit.i:             clamped := if(val <= -1, -0.99999, if(val >= 1, 0.99999, val));
./modules/tensor_mechanics/test/tests/rom_stress_update/lower_limit.i:             clamped := if(val <= -1, -0.99999, if(val >= 1, 0.99999, val));
./modules/tensor_mechanics/test/tests/rom_stress_update/lower_limit.i:             clamped := if(val <= -1, -0.99999, if(val >= 1, 0.99999, val));
./modules/tensor_mechanics/test/tests/finite_strain_elastic/elastic_rotation_test.i:    value = 'if(t<=1.0, delta*t, (1.0+delta)*cos(pi/2*(t-t0)) - 1.0)'
./modules/tensor_mechanics/test/tests/finite_strain_elastic/elastic_rotation_test.i:    value = 'if(t<=1.0, 0.0, (1.0+delta)*sin(pi/2*(t-t0)))'
./modules/tensor_mechanics/test/tests/finite_strain_elastic/elastic_rotation_test.i:    value = 'if(t<=1.0, delta*t, (1.0+delta)*cos(pi/2.0*(t-t0)) - sin(pi/2.0*(t-t0)) - 1.0)'
./modules/tensor_mechanics/test/tests/finite_strain_elastic/elastic_rotation_test.i:    value = 'if(t<=1.0, 0.0, cos(pi/2.0*(t-t0)) + (1+delta)*sin(pi/2.0*(t-t0)) - 1.0)'
./modules/tensor_mechanics/test/tests/finite_strain_elastic/elastic_rotation_test.i:    value = 'if(t<=1.0, 0.0, -sin(pi/2.0*(t-t0)))'
./modules/tensor_mechanics/test/tests/finite_strain_elastic/elastic_rotation_test.i:    value = 'if(t<=1.0, 0.0, cos(pi/2.0*(t-t0)) - 1.0)'
./modules/tensor_mechanics/test/tests/capped_weak_plane/small_deform6.i:    function = 'if(t<30,0.2*t,6)'
./modules/tensor_mechanics/test/tests/capped_weak_plane/small_deform6.i:    function = 'if(t<30,if(t<10,0,t),30-0.2*t)'
./modules/tensor_mechanics/test/tests/capped_weak_plane/small_deform6.i:    function = 'if(t<15,3*t,45)+if(t<30,0,45-3*t)'
./modules/tensor_mechanics/test/tests/capped_weak_plane/pull_push_h.i:    function = 'if(t>1,-2.0+t,-t)'
./modules/tensor_mechanics/test/tests/capped_weak_plane/pull_push.i:    function = 'if(t>1,-2.0+t,-t)'
./modules/tensor_mechanics/test/tests/weak_plane_tensile/large_deform2.i:    function = '0.70710678*y+0.70710678*z-y+if(t>0,1,0)'
./modules/tensor_mechanics/test/tests/weak_plane_tensile/large_deform2.i:    function = '-0.70710678*y+0.70710678*z-z+if(t>0,1,0)'
./modules/tensor_mechanics/test/tests/weak_plane_tensile/large_deform1.i:    function = '-1*y+0*z-z+if(t>0,0.5-y,0)' # note that this uses original nodal values of (x,y,z)
./modules/tensor_mechanics/test/tests/weak_plane_tensile/large_deform1.i:    function = '-1*y+0*z-z+if(t>0,0.5-y,0)' # note that this uses original nodal values of (x,y,z)
./modules/tensor_mechanics/test/tests/ad_smeared_cracking/cracking_rotation.i:    function = 'if(t<10,0,if(t>=100,1,1-cos((t-10)*pi/180)))'
./modules/tensor_mechanics/test/tests/ad_smeared_cracking/cracking_rotation.i:    function = 'if(t<10,0,if(t>=100,1,sin((t-10)*pi/180)))'
./modules/tensor_mechanics/test/tests/ad_smeared_cracking/cracking_rotation.i:    function = 'if(t<5,t*0.01,0.05-(t-5)*0.01)'
./modules/tensor_mechanics/test/tests/ad_smeared_cracking/cracking_function.i:    value = 'if(x > 0.667, 1.1e6, 1.2e6)'
./modules/tensor_mechanics/test/tests/ad_smeared_cracking/cracking_xyz.i:    value = 'if(t < 2, 0.00175, 0)'
./modules/tensor_mechanics/test/tests/capped_drucker_prager/random.i:    value = 'if(a<1E-3,0,a)'
./modules/tensor_mechanics/test/tests/capped_drucker_prager/random.i:    value = 'if(a<1E-3,0,a)'
./modules/tensor_mechanics/test/tests/capped_drucker_prager/random.i:    value = 'if(a<1E-3,0,a)'
./modules/tensor_mechanics/test/tests/recompute_radial_return/affine_plasticity.i:             if(t<=timeAtYield, -474*t,
./modules/tensor_mechanics/test/tests/recompute_radial_return/affine_plasticity.i:             if(t<=1, stressAtYield,
./modules/tensor_mechanics/test/tests/recompute_radial_return/affine_plasticity.i:             if(t<=timeAtYield, -474*t,
./modules/tensor_mechanics/test/tests/recompute_radial_return/affine_plasticity.i:             if(t<=1, stressAtYield,
./modules/tensor_mechanics/test/tests/recompute_radial_return/affine_plasticity.i:             if(t<=timeAtYield, 948*t,
./modules/tensor_mechanics/test/tests/recompute_radial_return/affine_plasticity.i:             if(t<=1, stressAtYield,
./modules/tensor_mechanics/test/tests/notched_plastic_block/cmc_planar.i:    function = 'if(plastic_strain>0,1,0)'
./modules/tensor_mechanics/test/tests/notched_plastic_block/biaxial_planar.i:    function = 'if(plastic_strain>0,1,0)'
./modules/tensor_mechanics/test/tests/interface_stress/multi.i:    value = 'r:=sqrt(x^2+y^2+z^2); if(r>1,0,1-3*r^2+2*r^3)'
./modules/tensor_mechanics/test/tests/interface_stress/multi.i:    value = 'r:=sqrt(x^2+y^2+z^2); 0.5-0.5*if(r>1,0,1-3*r^2+2*r^3)'
./modules/tensor_mechanics/test/tests/interface_stress/test.i:    value = 'r:=sqrt(x^2+y^2+z^2); if(r>1,0,1-3*r^2+2*r^3)'
./modules/tensor_mechanics/test/tests/mean_cap_TC/random02.i:    value = 'if(a<1E-3,0,a)'
./modules/tensor_mechanics/test/tests/mean_cap_TC/random03.i:    value = 'if(a<1E-3,0,a)'
./modules/tensor_mechanics/test/tests/mean_cap_TC/random01.i:    value = 'if(a<1E-3,0,a)'
./modules/tensor_mechanics/test/tests/mean_cap_TC/small_deform7.i:    function = 'if(t<1.5,-1E-7*x,1E-7*x)'
./modules/tensor_mechanics/test/tests/mean_cap_TC/small_deform7.i:    function = 'if(t<1.5,3E-7*y,1E-7*y)'
./modules/tensor_mechanics/test/tests/mean_cap_TC/small_deform7.i:    function = 'if(t<1.5,5E-7*z,4E-7*z)'
./modules/tensor_mechanics/test/tests/mean_cap_TC/random04.i:    value = 'if(a<1E-3,0,a)'
./modules/tensor_mechanics/test/tests/cohesive_zone_model/bilinear_mixed_scale_strength.i:    function = 'if(x<0.5,1,100)'
./modules/tensor_mechanics/test/tests/cohesive_zone_model/bilinear_mixed_scale_strength.i:    function = 'if(t<=0.3,t,if(t<=0.6,0.3-(t-0.3),0.6-t))'
./modules/tensor_mechanics/test/tests/cohesive_zone_model/bilinear_mixed.i:    function = 'if(t<=0.3,t,if(t<=0.6,0.3-(t-0.3),0.6-t))'
./modules/tensor_mechanics/test/tests/weak_plane_shear/large_deform4.i:    value = 'if(a<1E-3,0,a)'
./modules/tensor_mechanics/test/tests/weak_plane_shear/small_deform3.i:    value = 'if(a<1E-3,0,a)'
./modules/tensor_mechanics/test/tests/weak_plane_shear/small_deform2.i:    function = 'if(t<1E-6,0,3*t)'
./modules/tensor_mechanics/test/tests/weak_plane_shear/small_deform2.i:    function = 'if(t<1E-6,0,5*(t-0.01E-6))'
./modules/tensor_mechanics/test/tests/weak_plane_shear/small_deform2.i:    function = 'if(t<1E-6,t,2E-6-t)'
./modules/tensor_mechanics/test/tests/weak_plane_shear/small_deform_harden2.i:    function = 'if(t<1E-6,0,3*t)'
./modules/tensor_mechanics/test/tests/weak_plane_shear/small_deform_harden2.i:    function = 'if(t<1E-6,0,5*(t-0.01E-6))'
./modules/tensor_mechanics/test/tests/weak_plane_shear/small_deform_harden2.i:    function = 'if(t<1E-6,t,2E-6-t)'
./modules/tensor_mechanics/test/tests/weak_plane_shear/large_deform2.i:    function = '0.70710678*y+0.70710678*z-y+if(t>0,1,0)'
./modules/tensor_mechanics/test/tests/weak_plane_shear/large_deform2.i:    function = '-0.70710678*y+0.70710678*z-z+if(t>0,1,0)'
./modules/tensor_mechanics/test/tests/weak_plane_shear/large_deform3.i:    value = 'if(a<1E-3,0,a)'
./modules/tensor_mechanics/test/tests/weak_plane_shear/large_deform_harden3.i:    value = 'if(a<1E-3,0,a)'
./modules/tensor_mechanics/test/tests/weak_plane_shear/small_deform4.i:    function = 'if(t<1E-6,0,3*(t-1E-6)*(t-1E-6)*1E6)'
./modules/tensor_mechanics/test/tests/weak_plane_shear/small_deform4.i:    function = 'if(t<1E-6,0,5*(t-1E-6)*(t-1E-6)*1E6)'
./modules/tensor_mechanics/test/tests/weak_plane_shear/small_deform4.i:    function = 'if(t<1E-6,t,1E-6)'
./modules/tensor_mechanics/test/tests/weak_plane_shear/large_deform1.i:    function = '-1*y+0*z-z+if(t>0,0.5-y,0)' # note that this uses original nodal values of (x,y,z)
./modules/tensor_mechanics/test/tests/weak_plane_shear/large_deform1.i:    function = '-1*y+0*z-z+if(t>0,0.5-y,0)' # note that this uses original nodal values of (x,y,z)
./modules/tensor_mechanics/test/tests/mohr_coulomb/planar_hard4.i:    value = 'if((a<1E-5)&(b<1E-5)&(c<1E-5)&(d<1E-5)&(g<1E-5)&(h<1E-5),0,abs(a)+abs(b)+abs(c)+abs(d)+abs(g)+abs(h))'
./modules/tensor_mechanics/test/tests/mohr_coulomb/planar_hard4.i:    value = 'if(abs(a-b)<1E-6,0,1E6*abs(a-b))'
./modules/tensor_mechanics/test/tests/mohr_coulomb/planar_hard5.i:    value = 'if((a<1E-5)&(b<1E-5)&(c<1E-5)&(d<1E-5)&(g<1E-5)&(h<1E-5),0,abs(a)+abs(b)+abs(c)+abs(d)+abs(g)+abs(h))'
./modules/tensor_mechanics/test/tests/mohr_coulomb/planar_hard5.i:    value = 'if(abs(a-b)<1E-6,0,1E6*abs(a-b))'
./modules/tensor_mechanics/test/tests/mohr_coulomb/many_deforms_cap.i:    value = 'if(a<1E-3,0,a)'
./modules/tensor_mechanics/test/tests/mohr_coulomb/random_planar.i:    value = 'if(a<1E-3,0,a)'
./modules/tensor_mechanics/test/tests/mohr_coulomb/random.i:    value = 'if(a<1E-3,0,a)'
./modules/tensor_mechanics/test/tests/mohr_coulomb/planar_hard3.i:    value = 'if((a<1E-5)&(b<1E-5)&(c<1E-5)&(d<1E-5)&(g<1E-5)&(h<1E-5),0,abs(a)+abs(b)+abs(c)+abs(d)+abs(g)+abs(h))'
./modules/tensor_mechanics/test/tests/mohr_coulomb/planar_hard3.i:    value = 'if(abs(a-b)<1E-6,0,1E6*abs(a-b))'
./modules/tensor_mechanics/test/tests/mean_cap/random.i:    value = 'if(a<1E-3,0,a)'
./modules/tensor_mechanics/examples/coal_mining/cosserat_mc_wp_sticky.i:    value = 'if(y<ymin+(ymax-ymin)*min(t/end_t,1),minval,if(y<ymin+(ymax-ymin)*min(t/end_t,1)+slope,minval+(maxval-minval)*(y-(ymin+(ymax-ymin)*min(t/end_t,1)))/slope,maxval))'
./modules/tensor_mechanics/examples/coal_mining/cosserat_mc_wp_sticky.i:    value = 'if(y<ymin+(ymax-ymin)*min(t/end_t,1),minval,maxval)'
./modules/tensor_mechanics/examples/coal_mining/coarse.i:    value = 'if(y<ymin+(ymax-ymin)*min(t/end_t,1),minval,if(y<ymin+(ymax-ymin)*min(t/end_t,1)+slope,minval+(maxval-minval)*(y-(ymin+(ymax-ymin)*min(t/end_t,1)))/slope,maxval))'
./modules/tensor_mechanics/examples/coal_mining/coarse.i:    value = 'if(y<ymin+(ymax-ymin)*min(t/end_t,1),minval,maxval)'
./modules/tensor_mechanics/examples/coal_mining/cosserat_mc_wp_sticky_longitudinal.i:    value = 'if(y<ymin+(ymax-ymin)*min(t/end_t,1),minval,if(y<ymin+(ymax-ymin)*min(t/end_t,1)+slope,minval+(maxval-minval)*(y-(ymin+(ymax-ymin)*min(t/end_t,1)))/slope,maxval))'
./modules/tensor_mechanics/examples/coal_mining/cosserat_mc_wp_sticky_longitudinal.i:    value = 'if(y<ymin+(ymax-ymin)*min(t/end_t,1),minval,maxval)'
./modules/tensor_mechanics/examples/coal_mining/fine.i:    value = 'if(y<ymin+(ymax-ymin)*min(t/end_t,1),minval,if(y<ymin+(ymax-ymin)*min(t/end_t,1)+slope,minval+(maxval-minval)*(y-(ymin+(ymax-ymin)*min(t/end_t,1)))/slope,maxval))'
./modules/tensor_mechanics/examples/coal_mining/fine.i:    value = 'if(y<ymin+(ymax-ymin)*min(t/end_t,1),minval,maxval)'
./modules/richards/test/tests/user_objects/uo1.i:    value = if(x>ss,1,if(x<sn,0,kn+(((x-sn)/(ss-sn))^2)*(c-1)*(ks-kn)/(c-((x-sn)/(ss-sn)))))
./modules/richards/test/tests/user_objects/uo1.i:    value = if(x>ss,1,if(x<sn,0,kn+(((x-sn)/(ss-sn))^2)*(c-1)*(ks-kn)/(c-((x-sn)/(ss-sn)))))
./modules/richards/test/tests/user_objects/uo1.i:    value = if(x>ss,1,if(x<sn,0,kn+(((x-sn)/(ss-sn))^2)*(c-1)*(ks-kn)/(c-((x-sn)/(ss-sn)))))
./modules/richards/test/tests/user_objects/uo1.i:    value = if(x>simm,1,0)
./modules/richards/test/tests/user_objects/uo1.i:    value = if(x>simm,1,0)
./modules/richards/test/tests/user_objects/uo1.i:    value = if(x>simm,1,0)
./modules/richards/test/tests/user_objects/uo3.i:    value = if(x<pcut,scut+dscut*(x-pcut),(1+max((-x)*al,0)^(1/(1-m)))^(-m))
./modules/richards/test/tests/user_objects/uo3.i:    value = if(x<pcut,scut+dscut*(x-pcut),(1+max((-x)*al,0)^(1/(1-m)))^(-m))
./modules/richards/test/tests/user_objects/uo3.i:    value = if(x<pcut,scut+dscut*(x-pcut),(1+max((-x)*al,0)^(1/(1-m)))^(-m))
./modules/richards/test/tests/user_objects/uo2.i:    value = if(x>0,(0.00654576947608E-3*x+1.04357716547E-13*x^2),0)+if(x<0,0.1*(e^(6.54576947608E-5*x)-1),0)
./modules/richards/test/tests/user_objects/uo2.i:    value = if(x>0,(0.00654576947608E-3*x+1.04357716547E-13*x^2),0)+if(x<0,0.1*(e^(6.54576947608E-5*x)-1),0)
./modules/richards/test/tests/user_objects/uo2.i:    value = if(x>0,(0.00654576947608E-3*x+1.04357716547E-13*x^2),0)+if(x<0,0.1*(e^(6.54576947608E-5*x)-1),0)
./modules/richards/test/tests/user_objects/uo2.i:    value = if(x>0,-(molar_mass*(-2+(2*pow(2,0.3333333333333333)*(a-3*b*(b*x+rt)))/pow(-2*pow(a,3)+9*pow(a,2)*b*(-2*b*x+rt)+pow(pow(a,3)*(a*pow(2*a+9*b*(2*b*x-rt),2)-4*pow(a-3*b*(b*x+rt),3)),0.5),0.3333333333333333)+(pow(2,0.6666666666666666)*pow(-2*pow(a,3)+9*pow(a,2)*b*(-2*b*x+rt)+pow(pow(a,3)*(a*pow(2*a+9*b*(2*b*x-rt),2)-4*pow(a-3*b*(b*x+rt),3)),0.5),0.3333333333333333))/a))/(6.*b)+(molar_mass*(-2+(2*pow(2,0.3333333333333333)*(a-3*b*(b*0+rt)))/pow(-2*pow(a,3)+9*pow(a,2)*b*(-2*b*0+rt)+pow(pow(a,3)*(a*pow(2*a+9*b*(2*b*0-rt),2)-4*pow(a-3*b*(b*0+rt),3)),0.5),0.3333333333333333)+(pow(2,0.6666666666666666)*pow(-2*pow(a,3)+9*pow(a,2)*b*(-2*b*0+rt)+pow(pow(a,3)*(a*pow(2*a+9*b*(2*b*0-rt),2)-4*pow(a-3*b*(b*0+rt),3)),0.5),0.3333333333333333))/a))/(6.*b),infinityratio*molar_mass*(e^(slope0*x)-1))
./modules/richards/test/tests/user_objects/uo2.i:    value = if(x>0,-(molar_mass*(-2+(2*pow(2,0.3333333333333333)*(a-3*b*(b*x+rt)))/pow(-2*pow(a,3)+9*pow(a,2)*b*(-2*b*x+rt)+pow(pow(a,3)*(a*pow(2*a+9*b*(2*b*x-rt),2)-4*pow(a-3*b*(b*x+rt),3)),0.5),0.3333333333333333)+(pow(2,0.6666666666666666)*pow(-2*pow(a,3)+9*pow(a,2)*b*(-2*b*x+rt)+pow(pow(a,3)*(a*pow(2*a+9*b*(2*b*x-rt),2)-4*pow(a-3*b*(b*x+rt),3)),0.5),0.3333333333333333))/a))/(6.*b),infinityratio*molar_mass*(e^(slope0*x)-1))
./modules/richards/test/tests/user_objects/uo2.i:    value = if(x>0,-(molar_mass*(-2+(2*pow(2,0.3333333333333333)*(a-3*b*(b*x+rt)))/pow(-2*pow(a,3)+9*pow(a,2)*b*(-2*b*x+rt)+pow(pow(a,3)*(a*pow(2*a+9*b*(2*b*x-rt),2)-4*pow(a-3*b*(b*x+rt),3)),0.5),0.3333333333333333)+(pow(2,0.6666666666666666)*pow(-2*pow(a,3)+9*pow(a,2)*b*(-2*b*x+rt)+pow(pow(a,3)*(a*pow(2*a+9*b*(2*b*x-rt),2)-4*pow(a-3*b*(b*x+rt),3)),0.5),0.3333333333333333))/a))/(6.*b),infinityratio*molar_mass*(e^(slope0*x)-1))
./modules/richards/test/tests/user_objects/uo2.i:    value = if(x<zero_pt,0,if(x>cut_limit,dens0*exp(x/bulk_mod),(3*cut_limit-2*x-zero_pt)*(x-zero_pt)*(x-zero_pt)*dens0*exp(x/bulk_mod)/(cut_limit-zero_pt)/(cut_limit-zero_pt)/(cut_limit-zero_pt)))
./modules/richards/test/tests/user_objects/uo2.i:    value = if(x<zero_pt,0,if(x>cut_limit,dens0*exp(x/bulk_mod),(3*cut_limit-2*x-zero_pt)*(x-zero_pt)*(x-zero_pt)*dens0*exp(x/bulk_mod)/(cut_limit-zero_pt)/(cut_limit-zero_pt)/(cut_limit-zero_pt)))
./modules/richards/test/tests/user_objects/uo2.i:    value = if(x<zero_pt,0,if(x>cut_limit,dens0*exp(x/bulk_mod),(3*cut_limit-2*x-zero_pt)*(x-zero_pt)*(x-zero_pt)*dens0*exp(x/bulk_mod)/(cut_limit-zero_pt)/(cut_limit-zero_pt)/(cut_limit-zero_pt)))
./modules/richards/test/tests/buckley_leverett/bl20_lumped.i:    value = 1000000*(1-min(x/5,1))-if(x<5,0,300000)
./modules/richards/test/tests/buckley_leverett/bl20_lumped_fu.i:    value = 1000000*(1-min(x/5,1))-if(x<5,0,300000)
./modules/richards/test/tests/buckley_leverett/bl22_lumped.i:    value = 1000000*(1-min(x/5,1))-if(x<5,0,100000)
./modules/richards/test/tests/buckley_leverett/bl22_lumped_fu.i:    value = 1000000*(1-min(x/5,1))-if(x<5,0,100000)
./modules/fluid_properties/test/tests/co2/co2.i:    value = if(x<1,280,if(x<2,360,500))
./modules/fluid_properties/test/tests/brine/brine_tabulated.i:    value = 'if(x<2,20e6, 40e6)'
./modules/fluid_properties/test/tests/brine/brine_tabulated.i:    value = 'if(x<1, 323.15, 473.15)'
./modules/fluid_properties/test/tests/brine/brine_tabulated.i:    value = 'if(x<2,0.1047, 0.2261)'
./modules/fluid_properties/test/tests/brine/brine.i:    value = 'if(x<2,20e6, 40e6)'
./modules/fluid_properties/test/tests/brine/brine.i:    value = 'if(x<1, 323.15, 473.15)'
./modules/fluid_properties/test/tests/brine/brine.i:    value = 'if(x<2,0.1047, 0.2261)'
./modules/fluid_properties/test/tests/water/water.i:    value = 'if(x<2, 300, 500)'
./modules/fluid_properties/test/tests/water/water.i:    value = 'if(x<1,3e6, if(x<2, 80e6, 3e6))'
./modules/porous_flow/test/tests/flux_limited_TVD_advection/fltvd_2D_blocks.i:    function = 'if(x<0.1 | x>0.3, 0, 1)'
./modules/porous_flow/test/tests/flux_limited_TVD_advection/fltvd_2D_blocks.i:    function = 'if(x<0.7 | x > 0.9, 0, 1)'
./modules/porous_flow/test/tests/flux_limited_TVD_advection/fltvd_2D.i:    function = 'if(x<0.1,0,if(x>0.3,0,1))'
./modules/porous_flow/test/tests/flux_limited_TVD_advection/fltvd_3D.i:    function = 'if(x<0.1,0,if(x>0.3,0,1))'
./modules/porous_flow/test/tests/flux_limited_TVD_advection/fltvd_1D_adaptivity.i:    function = 'if(x<0.1,0,if(x>0.3,0,1))'
./modules/porous_flow/test/tests/flux_limited_TVD_advection/fltvd_1D.i:    function = 'if(x<0.1,0,if(x>0.3,0,1))'
./modules/porous_flow/test/tests/flux_limited_TVD_advection/fltvd_2D_angle.i:    function = 'if(x<0.1 | x > 0.3 | y < 0.1 | y > 0.3, 0, 1)'
./modules/porous_flow/test/tests/flux_limited_TVD_advection/fltvd_2D_trimesh.i:    function = 'if(x<0.1,0,if(x>0.305,0,1))'
./modules/porous_flow/test/tests/recover/pffltvd.i:    function = 'if(x<0.1,0,if(x>0.3,0,1))'
./modules/porous_flow/test/tests/hysteresis/2phasePS_relperm_2.i:    function = 'if(t <= 15, 20, -20)'
./modules/porous_flow/test/tests/hysteresis/hys_order_02.i:    value = '30 * if(t <= 7, -1, if(t <= 10, 1, if(t <= 12, -1, 1)))'
./modules/porous_flow/test/tests/hysteresis/hys_order_05.i:    value = '30 * if(t <= 2, -1, if(t <= 7, 1, -1))'
./modules/porous_flow/test/tests/hysteresis/1phase_relperm.i:    function = 'if(t <= 5, -10, 10)'
./modules/porous_flow/test/tests/hysteresis/hys_order_07.i:    value = '30 * if(t <= 1, -1, 1)'
./modules/porous_flow/test/tests/hysteresis/2phasePP.i:    function = 'if(t <= 9, 10, -10)'
./modules/porous_flow/test/tests/hysteresis/hys_order_01.i:    value = '30 * if(t <= 4, -1, if(t <= 7, 1, -1))'
./modules/porous_flow/test/tests/hysteresis/1phase_relperm_2.i:    function = 'if(t <= 3, -10, if(t <= 5, 10, if(t <= 13, -10, 10)))'
./modules/porous_flow/test/tests/hysteresis/hys_order_03.i:    value = '30 * if(t <= 8, -1, if(t <= 15, 1, if(t <= 20, -1, if(t <= 24, 1, if(t <= 27, -1, if(t <= 30, 1, -1))))))'
./modules/porous_flow/test/tests/hysteresis/2phasePS.i:    function = 'if(t <= 9, 10, -10)'
./modules/porous_flow/test/tests/hysteresis/1phase.i:    function = 'if(t <= 9, -10, 10)'
./modules/porous_flow/test/tests/hysteresis/hys_order_09.i:  value = '30 * if(t <= 1, -2, if(t <= 2, 1.5, -1))'
./modules/porous_flow/test/tests/hysteresis/2phasePS_relperm.i:    function = 'if(t <= 9, 10, -10)'
./modules/porous_flow/test/tests/hysteresis/2phasePP_2.i:  function = 'if(t <= 14, 10, if(t <= 25, -10, 10))'
./modules/porous_flow/test/tests/hysteresis/1phase_3rd.i:    function = 'if(t <= 9, -10, if(t <= 16, 10, if(t <= 22, -10, 10)))'
./modules/porous_flow/test/tests/hysteresis/hys_order_08.i:  value = '30 * if(t <= 1, -2, if(t <= 2, 1.5, -2))'
./modules/porous_flow/test/tests/hysteresis/2phasePS_2.i:  function = 'if(t <= 14, 10, if(t <= 25, -10, 10))'
./modules/porous_flow/test/tests/dirackernels/pls03_action.i:    #function = if((x<1)&(y<0.5),1E7,-1E7)
./modules/porous_flow/test/tests/dirackernels/pls03_action.i:    function = if((x<1)&(y>0.5),1E7,-1E7)
./modules/porous_flow/test/tests/dirackernels/pls03_action.i:    #function = if((x>1)&(y<0.5),1E7,-1E7)
./modules/porous_flow/test/tests/dirackernels/pls03_action.i:    #function = if((x>1)&(y>0.5),1E7,-1E7)
./modules/porous_flow/test/tests/dirackernels/pls03.i:    #function = if((x<1)&(y<0.5),1E7,-1E7)
./modules/porous_flow/test/tests/dirackernels/pls03.i:    function = if((x<1)&(y>0.5),1E7,-1E7)
./modules/porous_flow/test/tests/dirackernels/pls03.i:    #function = if((x>1)&(y<0.5),1E7,-1E7)
./modules/porous_flow/test/tests/dirackernels/pls03.i:    #function = if((x>1)&(y>0.5),1E7,-1E7)
./modules/porous_flow/test/tests/numerical_diffusion/fltvd.i:    function = 'if(x<0.1,0,if(x>0.3,0,1))'
./modules/porous_flow/test/tests/numerical_diffusion/fltvd_no_antidiffusion.i:    function = 'if(x<0.1,0,if(x>0.3,0,1))'
./modules/porous_flow/test/tests/numerical_diffusion/pffltvd_action.i:    function = 'if(x<0.1,0,if(x>0.3,0,1))'
./modules/porous_flow/test/tests/numerical_diffusion/pffltvd.i:    function = 'if(x<0.1,0,if(x>0.3,0,1))'
./modules/porous_flow/test/tests/numerical_diffusion/fully_saturated_action.i:    function = 'if(x<0.1,0,if(x>0.3,0,1))'
./modules/porous_flow/test/tests/numerical_diffusion/framework.i:    function = 'if(x<0.1,0,if(x>0.3,0,1))'
./modules/porous_flow/test/tests/numerical_diffusion/fltvd_none.i:    function = 'if(x<0.1,0,if(x>0.3,0,1))'
./modules/porous_flow/test/tests/numerical_diffusion/no_action.i:    function = 'if(x<0.1,0,if(x>0.3,0,1))'
./modules/porous_flow/test/tests/sinks/s06.i:    value = 'vol*por*dens0*exp(pp/bulk)*if(pp>=0,1,pow(1+pow(-al*pp,1.0/(1-m)),-m))'
./modules/porous_flow/test/tests/sinks/s06.i:    value = 'fcn*if(pp>center,m,if(pp<themin,0,m/c/c/c*(2*(pp-center)+c)*((pp-center)-c)*((pp-center)-c)))'
./modules/porous_flow/test/tests/sinks/s06.i:    value = 'vol*por*dens0*exp(pp/bulk)*if(pp>=0,1,pow(1+pow(-al*pp,1.0/(1-m)),-m))'
./modules/porous_flow/test/tests/sinks/s06.i:    value = 'fcn*if(pp>center,m,if(pp<themin,0,m/c/c/c*(2*(pp-center)+c)*((pp-center)-c)*((pp-center)-c)))'
./modules/porous_flow/test/tests/sinks/s04.i:    value = 'fcn*if(pp>0.8,1,if(pp<0.3,0.5,0.2+pp))'
./modules/porous_flow/test/tests/sinks/s04.i:    value = 'fcn*if(pp>0.8,1,if(pp<0.3,0.5,0.2+pp))'
./modules/porous_flow/test/tests/sinks/s05.i:    value = 'vol*por*dens0*exp(pp/bulk)*if(pp>=0,1,pow(1+pow(-al*pp,1.0/(1-m)),-m))'
./modules/porous_flow/test/tests/sinks/s05.i:    value = 'if(pp>center,fcn,fcn*exp(-0.5*(pp-center)*(pp-center)/sd/sd))'
./modules/porous_flow/test/tests/sinks/s05.i:    value = 'vol*por*dens0*exp(pp/bulk)*if(pp>=0,1,pow(1+pow(-al*pp,1.0/(1-m)),-m))'
./modules/porous_flow/test/tests/sinks/s05.i:    value = 'if(pp>center,fcn,fcn*exp(-0.5*(pp-center)*(pp-center)/sd/sd))'
./modules/porous_flow/test/tests/flux_limited_TVD_pflow/pffltvd_1D_adaptivity.i:    function = 'if(x<0.1,0,if(x>0.3,0,1))'
./modules/porous_flow/test/tests/flux_limited_TVD_pflow/pffltvd_1D.i:    function = 'if(x<0.1,0,if(x>0.3,0,1))'
./modules/porous_flow/test/tests/flux_limited_TVD_pflow/pffltvd_2D.i:    function = 'if(x<0.1,0,if(x>0.3,0,1))'
./modules/porous_flow/test/tests/flux_limited_TVD_pflow/pffltvd_2D_angle.i:    function = 'if(x<0.1 | x > 0.3 | y < 0.1 | y > 0.3, 0, 1)'
./modules/porous_flow/test/tests/flux_limited_TVD_pflow/pffltvd_2D_trimesh.i:    function = 'if(x<0.1,0,if(x>0.305,0,1))'
./modules/porous_flow/test/tests/flux_limited_TVD_pflow/jacobian_05.i:    function = 'if(x<1,0,if(x<4,sin(x-1),1))'
./modules/porous_flow/test/tests/flux_limited_TVD_pflow/pffltvd_3D.i:    function = 'if(x<0.1,0,if(x>0.3,0,1))'
./modules/porous_flow/test/tests/basic_advection/1phase.i:    function = 'if(x<0.1,1,0)'
./modules/porous_flow/test/tests/basic_advection/except1.i:    function = 'if(x<0.1,1,0)'
./modules/porous_flow/test/tests/basic_advection/2phase.i:    function = 'if(x<0.1,1,0)'
./modules/porous_flow/test/tests/basic_advection/except2.i:    function = 'if(x<0.1,1,0)'
./modules/porous_flow/test/tests/poro_elasticity/terzaghi_basicthm.i:    function = if(0.5*t<0.1,0.5*t,0.1)
./modules/porous_flow/test/tests/poro_elasticity/mandel_fully_saturated.i:    function = if(0.15*t<0.01,0.15*t,0.01)
./modules/porous_flow/test/tests/poro_elasticity/terzaghi_fully_saturated_volume.i:    function = if(0.5*t<0.1,0.5*t,0.1)
./modules/porous_flow/test/tests/poro_elasticity/terzaghi_constM.i:    function = if(0.5*t<0.1,0.5*t,0.1)
./modules/porous_flow/test/tests/poro_elasticity/mandel_basicthm.i:    function = if(0.15*t<0.01,0.15*t,0.01)
./modules/porous_flow/test/tests/poro_elasticity/terzaghi.i:    function = if(0.5*t<0.1,0.5*t,0.1)
./modules/porous_flow/test/tests/poro_elasticity/mandel.i:    function = if(0.15*t<0.01,0.15*t,0.01)
./modules/porous_flow/test/tests/poro_elasticity/mandel_constM.i:    function = if(0.15*t<0.01,0.15*t,0.01)
./modules/porous_flow/test/tests/poro_elasticity/mandel_fully_saturated_volume.i:    function = if(0.15*t<0.01,0.15*t,0.01)
./modules/porous_flow/examples/multiapp_fracture_flow/diffusion_multiapp/fracture_app_heat.i:    function = 'if(x<1E-6, 2, 0)'  # delta function
./modules/porous_flow/examples/multiapp_fracture_flow/diffusion_multiapp/fracture_app.i:    function = 'if(x<1E-6, 2, 0)'  # delta function
./modules/porous_flow/examples/multiapp_fracture_flow/diffusion_multiapp/single_var.i:    function = 'if(x<0.5, 2, 0)'  # delta function
./modules/porous_flow/examples/multiapp_fracture_flow/diffusion_multiapp/two_vars.i:    function = 'if(x<0.5, 2, 0)'  # delta function
./modules/porous_flow/examples/multiapp_fracture_flow/3dFracture/fracture_only_aperture_changing.i:    function = 'if(enclosing_element_normal_length = 0, 0, h_s * enclosing_element_normal_thermal_cond * 2 * enclosing_element_normal_length / (h_s * enclosing_element_normal_length * enclosing_element_normal_length + enclosing_element_normal_thermal_cond * 2 * enclosing_element_normal_length))'
./modules/porous_flow/examples/tutorial/06_KT.i:    function = '0.5*if(x*x+y*y<1.01,1,0)'
./modules/porous_flow/examples/tutorial/06.i:    function = '0.5*if(x*x+y*y<1.01,1,0)'
./modules/porous_flow/examples/tidal/atm_tides.i:    value = 'if(t>0,5000 * sin(2 * pi * t / 3600.0 / 24.0),0)'
./modules/porous_flow/examples/tidal/atm_tides.i:    value = '-if(t>0,5000 * sin(2 * pi * t / 3600.0 / 24.0),0)'
./modules/porous_flow/examples/tidal/atm_tides_open_hole.i:    value = 'if(t>0,5000 * sin(2 * pi * t / 3600.0 / 24.0),0)'
./modules/porous_flow/examples/tidal/atm_tides_open_hole.i:    value = '-10000*z + if(t>0,5000 * sin(2 * pi * t / 3600.0 / 24.0),0)'
./modules/porous_flow/examples/tidal/atm_tides_open_hole.i:    value = '-if(t>0,5000 * sin(2 * pi * t / 3600.0 / 24.0),0)'
./modules/porous_flow/examples/coal_mining/coarse_with_fluid.i:    value = 'if(y<ymin+(ymax-ymin)*min(t/end_t,1),minval,if(y<ymin+(ymax-ymin)*min(t/end_t,1)+slope,minval+(maxval-minval)*(y-(ymin+(ymax-ymin)*min(t/end_t,1)))/slope,maxval))'
./modules/porous_flow/examples/coal_mining/coarse_with_fluid.i:    value = 'if(y<ymin+(ymax-ymin)*min(t/end_t,1),minval,maxval)'
./modules/porous_flow/examples/coal_mining/coarse_with_fluid.i:    value = 'if(y<ymin+(ymax-ymin)*min(t/end_t,1),maxval,minval)'
./modules/porous_flow/examples/coal_mining/fine_with_fluid.i:    value = 'if(y<ymin+(ymax-ymin)*min(t/end_t,1),minval,if(y<ymin+(ymax-ymin)*min(t/end_t,1)+slope,minval+(maxval-minval)*(y-(ymin+(ymax-ymin)*min(t/end_t,1)))/slope,maxval))'
./modules/porous_flow/examples/coal_mining/fine_with_fluid.i:    value = 'if(y<ymin+(ymax-ymin)*min(t/end_t,1),minval,maxval)'
./modules/porous_flow/examples/coal_mining/fine_with_fluid.i:    value = 'if(y<ymin+(ymax-ymin)*min(t/end_t,1),maxval,minval)'
./modules/porous_flow/examples/ates/ates.i:    value = 'if(t >= ${start_injection1} & t < ${end_injection1}, 1,
./modules/porous_flow/examples/ates/ates.i:             if(t >= ${start_injection2} & t < ${end_injection2}, 1,
./modules/porous_flow/examples/ates/ates.i:             if(t >= ${start_injection3} & t < ${end_injection3}, 1,
./modules/porous_flow/examples/ates/ates.i:             if(t >= ${start_injection4} & t < ${end_injection4}, 1,
./modules/porous_flow/examples/ates/ates.i:             if(t >= ${start_injection5} & t < ${end_injection5}, 1,
./modules/porous_flow/examples/ates/ates.i:             if(t >= ${start_injection6} & t < ${end_injection6}, 1,
./modules/porous_flow/examples/ates/ates.i:             if(t >= ${start_injection7} & t < ${end_injection7}, 1,
./modules/porous_flow/examples/ates/ates.i:             if(t >= ${start_injection8} & t < ${end_injection8}, 1,
./modules/porous_flow/examples/ates/ates.i:             if(t >= ${start_injection9} & t < ${end_injection9}, 1,
./modules/porous_flow/examples/ates/ates.i:             if(t >= ${start_injection10} & t < ${end_injection10}, 1, 0))))))))))'
./modules/porous_flow/examples/ates/ates.i:    value = 'if(t >= ${start_production1} & t < ${end_production1}, 1,
./modules/porous_flow/examples/ates/ates.i:             if(t >= ${start_production2} & t < ${end_production2}, 1,
./modules/porous_flow/examples/ates/ates.i:             if(t >= ${start_production3} & t < ${end_production3}, 1,
./modules/porous_flow/examples/ates/ates.i:             if(t >= ${start_production4} & t < ${end_production4}, 1,
./modules/porous_flow/examples/ates/ates.i:             if(t >= ${start_production5} & t < ${end_production5}, 1,
./modules/porous_flow/examples/ates/ates.i:             if(t >= ${start_production6} & t < ${end_production6}, 1,
./modules/porous_flow/examples/ates/ates.i:             if(t >= ${start_production7} & t < ${end_production7}, 1,
./modules/porous_flow/examples/ates/ates.i:             if(t >= ${start_production8} & t < ${end_production8}, 1,
./modules/porous_flow/examples/ates/ates.i:             if(t >= ${start_production9} & t < ${end_production9}, 1,
./modules/porous_flow/examples/ates/ates.i:             if(t >= ${start_production10} & t < ${end_production10}, 1, 0))))))))))'
./modules/porous_flow/examples/co2_intercomparison/1Dradial/properties.i:    value = 'if(x<1,12e6,if(x<2,16e6,if(x<3,20e6,24e6)))'
./modules/porous_flow/examples/groundwater/ex01.i:    value = 'if(z <= -90, low, if(z >= -80, up, (up * (z + 90) - low * (z + 80)) / (10.0)))'
