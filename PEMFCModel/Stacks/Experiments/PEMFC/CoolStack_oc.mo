within PEMFCModel.Stacks.Experiments.PEMFC;

model CoolStack_oc "Cooled PEMFC stack test model"
    
  extends .Modelon.Icons.Experiment;

  replaceable package Medium_an = .PEMFCModel4par.Medium.NASASteamH2N2O2 constrainedby
    .FuelCell.Media.Templates.ReactionGas "Anodic medium";
  replaceable package Medium_cath = .FuelCell.Media.PreDefined.IdealGases.NASAMoistAir constrainedby
    .FuelCell.Media.Templates.ReactionGas "Cathodic medium";
  package MediumWater = .Modelon.Media.PreDefined.TwoPhase.WaterIF97 "Water medium";
  
  parameter .Modelica.Units.SI.MassFlowRate m_flow_an_min = 2.6e-4 "Minimum value of anodic gas inlet mass flow rate" annotation(Dialog(group = "Anode side"));
  parameter .Modelica.Units.SI.MassFlowRate m_flow_an_max = m_flow_an_min "Maximum value of anodic gas inlet mass flow rate" annotation(Dialog(group = "Anode side"));
  parameter .Modelica.Units.SI.Pressure p_an_in = 1.4e5 "Start inlet pressure at anode side" annotation(Dialog(group = "Anode side"));
  parameter .Modelica.Units.SI.Pressure p_an_out = 1.3e5 "Outlet pressure at anode side" annotation(Dialog(group = "Anode side"));
  parameter .Modelica.Units.SI.Temperature T_an_in = 70+273.15 "Anodic gas inlet temperature" annotation(Dialog(group = "Anode side"));
  parameter .Modelica.Units.SI.MassFraction X_an[Medium_an.nS] = Medium_an.moleToMassFractions({0.629,0.101,0.27,0},Medium_an.MMX) "Anodic gas inlet mass-based composition"  annotation(Dialog(group = "Anode side"));
  
  parameter .Modelica.Units.SI.MassFlowRate m_flow_cath_min = 1.376e-3 "Minimum value of cathodic gas inlet mass flow rate" annotation(Dialog(group = "Cathode side"));
  parameter .Modelica.Units.SI.MassFlowRate m_flow_cath_max = m_flow_cath_min "Maximum value of cathodic gas inlet mass flow rate" annotation(Dialog(group = "Cathode side"));
  parameter .Modelica.Units.SI.Pressure p_cath_in = 1.2e5 "Start inlet pressure at cathode side" annotation(Dialog(group = "Cathode side"));
  parameter .Modelica.Units.SI.Pressure p_cath_out = 1.1e5 "Outlet pressure at cathode side" annotation(Dialog(group = "Cathode side"));
  parameter .Modelica.Units.SI.Temperature T_cath_in = 70+273.15 "Cathodic gas inlet temperature" annotation(Dialog(group = "Cathode side"));
  parameter .Modelica.Units.SI.MassFraction X_cath[Medium_cath.nS] = Medium_cath.moleToMassFractions({0,0,0.117,0.185,0.697},Medium_cath.MMX) "Cathodic gas inlet mass-based composition" annotation(Dialog(group = "Cathode side"));

  parameter .Modelica.Units.SI.MassFlowRate m_flow_water_min  = 1e-3 "Minimum value of water inlet mass flow rate" annotation(Dialog(group = " Cooling"));
  parameter .Modelica.Units.SI.MassFlowRate m_flow_water_max = m_flow_water_min "Maximum value of water inlet mass flow rate" annotation(Dialog(group = " Cooling"));
  parameter .Modelica.Units.SI.Temperature T_water_in = 65+273.15 "Water inlet temperature" annotation(Dialog(group = " Cooling"));
  parameter .Modelica.Units.SI.Pressure p_water_in = 1.6e5 "Start water inlet pressure" annotation(Dialog(group = " Cooling"));
  parameter .Modelica.Units.SI.Pressure p_water_out = 1.2e5 "Water outlet pressure" annotation(Dialog(group = " Cooling"));
  //parameter .Modelica.Units.SI.MassFlowRate m_out=m_air*X_air[cath[1]]*(1 - 1/cath_stoich)
    //"calculated mass flow rate of oxygen gas at cathode outlet, may be used for control of the inlet cathode flow rate";
  parameter .Modelica.Units.SI.Current I_min = 2 "Minimum value of load current" annotation(Dialog(group = "Load"));
  parameter .Modelica.Units.SI.Current I_max = 225 "Maximum value of load current" annotation(Dialog(group = "Load", enable = not constant_current));
  parameter Boolean constant_current = true "Opportunely settings of ramp input signals" annotation(Dialog(group = "Load"));
 
  parameter .Modelica.Units.SI.Time simulation_time = 3600 "Simulation time" annotation(Dialog(tab = "Simulation options"));  
 
  .FuelCell.Sources.GasFlowBoundary flowAnode(
    redeclare package Medium = Medium_an,
    X = X_an,
    T = T_an_in,
    use_flow_in = true,
    use_Th_in = false) annotation (Placement(transformation(
          extent={{-101.93906632502666,19.768425703538625},{-81.93906632502666,39.768425703538625}},
                                     rotation=0.0,origin = {0.0,0.0})));

  .FuelCell.Sources.GasPressureBoundary sinkAnode(
    X = X_an,
    p = p_an_out,
    T = T_an_in,
    redeclare package Medium = Medium_an) annotation (Placement(transformation(extent={{99.22993225215629,19.6149661260782},{79.22993225215629,39.6149661260782}},
          rotation=0.0,origin = {0.0,0.0})));

  .Modelica.Electrical.Analog.Basic.Ground ground annotation (Placement(transformation(extent={{19.040856259459844,59.87114738360049},{31.040856259459844,71.87114738360049}},   rotation=
           0.0,origin = {0.0,0.0})));
        
  .Modelica.Electrical.Analog.Sources.RampCurrent current(
    startTime = if constant_current then simulation_time else 0,
    duration = if constant_current then 0 else simulation_time,
    offset = I_min,
    I = if constant_current then 0 else I_max - I_min) annotation (Placement(transformation(
        origin={0.0,64.0},
        extent={{7.0,-7.0},{-7.0,7.0}},
        rotation=-180.0)));
        
  .FuelCell.Sources.GasFlowBoundary flowCathode(
    redeclare package Medium = Medium_cath,
    X = X_cath,
    T = T_cath_in,
    use_flow_in = true,
    use_Th_in = false) annotation (Placement(transformation(
          extent={{-102.0,-30.0},{-82.0,-10.0}},
                                      rotation=0.0,origin = {0.0,0.0})));
        
  .FuelCell.Sources.GasPressureBoundary sinkCathode(
    redeclare package Medium = Medium_cath,
    p = p_cath_out,
    T = T_cath_in,X = X_cath) annotation (Placement(transformation(extent={{100.0,-30.0},{80.0,-10.0}},
          rotation=0.0,origin = {0.0,0.0})));
        
  .PEMFCModel.Stacks.PEMFC.Empirical.CoolStack coolStack(
    redeclare package Medium_an = Medium_an,
    redeclare package Medium_cath = Medium_cath,
    redeclare package Medium_cooling = MediumWater,
    X_feed_an = X_an,
    X_feed_cath = X_cath, X_start_anode = flowAnode.X,m_flow_start_anode = m_flow_an_min,p_start_out_anode = sinkAnode.p,p_start_in_anode = p_an_in,X_start_cathode = flowCathode.X,m_flow_start_cathode = m_flow_cath_min,p_start_out_cathode = sinkCathode.p,p_start_in_cathode = p_cath_in,T_start_in_anode = flowAnode.T,T_start_in_cathode = flowCathode.T,T_start_out_cathode = sinkCathode.T,T_start_out_anode = sinkAnode.T,positiveFlow_cooling = true,from_dp_cooling = false,positiveFlow_anode = true,positiveFlow_cathode = true,from_dp_anode = false,from_dp_cathode = false,p_start_in_cooling = p_water_in,T_start_in_cooling = sourceW.T,T_start_out_cooling = sinkP.T,m_flow_start_cooling = m_flow_water_min, h_inflow_an = flowAnode.fluidPort.h_outflow, h_inflow_cath = flowCathode.fluidPort.h_outflow, h_inflow_cooling = sourceW.fluidPort.h_outflow,p_start_out_cooling = sinkP.p,enable_setting = true,dp_smooth = 1e-2,mflow_smooth = 1e-8,dp_smooth_cooling = 1e-2,mflow_smooth_cooling = 1e-8,redeclare replaceable model Friction_anode = .Modelon.ThermoFluid.FlowChannels.PipeResistances.SinglePhase.LinearOperatingPointLoss(d0 = 0.5,dp0 = 0.1e5,m_flow0 = m_flow_an_min),redeclare replaceable model Friction_cathode = .Modelon.ThermoFluid.FlowChannels.PipeResistances.SinglePhase.LinearOperatingPointLoss(d0 = 1.2,dp0 = 0.1e5,m_flow0 = m_flow_cath_min),redeclare replaceable model HeatTransfer_anode = .Modelon.ThermoFluid.FlowChannels.HeatTransfer.SinglePhase.SinglePhase,redeclare replaceable model HeatTransfer_cathode = .Modelon.ThermoFluid.FlowChannels.HeatTransfer.SinglePhase.SinglePhase,redeclare replaceable model Friction_cooling = .Modelon.ThermoFluid.FlowChannels.PipeResistances.TwoPhase.DensityProfileFriction(p0_in = p_water_in,p0_out = p_water_out,mflow0 = m_flow_water_min),redeclare replaceable model HeatTransfer_cooling = .Modelon.ThermoFluid.FlowChannels.HeatTransfer.TwoPhase.SimpleTwoPhase,N = 5,A_cell = 300e-4,L_anode = 0.0813,D_anode = 0.023,D_cathode = 0.0445,L_cathode = 0.0813,L_cooling = 0.0813,D_cooling = 0.0375,diffusiveSpecies = {"O2","N2"},n_cell = 11,M_stack = 10.1,z = 75e-6,m_conc = 3.4e-4,alpha = 0.5,j_loss = 100,j_0 = 10,c1 = 5.75e-3) annotation (Placement(transformation(extent={{-22.0,0.0},{18.0,40.0}},rotation = 0.0,origin = {0.0,0.0})));

   .FuelCell.Sources.WaterFlowBoundary sourceW(redeclare package Medium = MediumWater, T = T_water_in,use_flow_in = true,use_Th_in = false) annotation (Placement(transformation(extent={{-58.0,-50.0},{-38.0,-30.0}},rotation = 0.0,origin = {0.0,0.0})));
        
  .FuelCell.Sources.WaterPressureBoundary sinkP(redeclare package Medium = MediumWater, T = T_water_in,p = p_water_out) annotation (Placement(transformation(extent={{-10.0,10.0},{10.0,-10.0}},
        rotation=180.0,
        origin={48.0,-40.0})));
        
  .Modelon.Visualizers.RealValue display_V(number = current.p.v - current.n.v, precision = 3) annotation (Placement(transformation(extent={{-90,-86},{-70,-66}})));
  .Modelon.Visualizers.RealValue display_P(number = current.p.i*(current.p.v - current.n.v), precision = 0) annotation (Placement(transformation(extent={{-66,-86},{-46,-66}})));
  .Modelon.Visualizers.RealValue display_I(number = current.p.i, precision = 2) annotation (Placement(transformation(extent={{-90,-102},{-70,-82}})));
  .Modelon.Visualizers.RealValue display_j(number = coolStack.summary.j_external*1e-4*1e3, precision = 1) annotation (Placement(transformation(extent={{-38,-86},{-18,-66}})));
  .Modelon.Visualizers.RealValue display_Q(number = sum(coolStack.coolingPipe.q.Q_flow), precision = 1) annotation (Placement(transformation(extent={{-38,-102},{-18,-82}})));
  .Modelon.Visualizers.RealValue display_T(precision = 1, number = .Modelica.Units.Conversions.to_degC(coolStack.subStack.cell.T_cell_avg)) annotation (Placement(transformation(extent={{-66,-102},{-46,-82}})));
  
  .FuelCell.Sensors.ReformateLongCompositionDisplay display_an_in(redeclare package Medium = Medium_an) annotation (Placement(transformation(extent={{-52.462753249644805,42.2554393386981},{-28.462753249644805,62.2554393386981}},rotation = 0.0,origin = {0.0,0.0})));
  .FuelCell.Sensors.ReformateLongCompositionDisplay display_an_out(redeclare package Medium = Medium_an) annotation (Placement(transformation(extent={{50,42},{74,62}})));
  .FuelCell.Sensors.AirCompositionDisplay display_cath_in(redeclare package Medium = Medium_cath) annotation (Placement(transformation(extent={{-53.744560661301904,-14.0},{-31.744560661301904,4.0}},rotation = 0.0,origin = {0.0,0.0})));
  .FuelCell.Sensors.AirCompositionDisplay display_cath_out(redeclare package Medium = Medium_cath) annotation (Placement(transformation(extent={{51.00000000000004,-13.0},{73.00000000000004,5.0}},rotation = 0.0,origin = {0.0,0.0})));    
  
  .FuelCell.Sensors.MultiDisplay_phTmdot display_phTmdot_an_in(displayUnits = true, redeclare package Medium = Medium_an, precision_p = 5, precision_mdot = 2) annotation (Placement(transformation(extent={{-78,30},{-58,50}})));
  .FuelCell.Sensors.MultiDisplay_phTmdot display_phTmdot_an_out(displayUnits = true, redeclare package Medium = Medium_an) annotation (Placement(transformation(extent={{24,30},{44,50}})));
  .FuelCell.Sensors.MultiDisplay_phTmdot display_phTmdot_cath_in(displayUnits = true, redeclare package Medium = Medium_cath,precision_p = 5, precision_mdot = 2) annotation (Placement(transformation(extent={{-78.0,-20.0},{-58.0,0.0}},rotation = 0.0,origin = {0.0,0.0})));
  .FuelCell.Sensors.MultiDisplay_phTmdot display_phTmdot_cath_out(displayUnits = true, redeclare package Medium = Medium_cath) annotation (Placement(transformation(extent={{26,-20},{46,0}})));
  .FuelCell.Sensors.MultiDisplay_phTmdot display_phTmdot_w_in(displayUnits = true, redeclare package Medium = MediumWater) annotation (Placement(transformation(extent={{-35.547635814819174,-42.0},{-15.547635814819174,-22.0}},rotation = 0.0,origin = {0.0,0.0})));
  .FuelCell.Sensors.MultiDisplay_phTmdot display_phTmdot_w_out(displayUnits = true, redeclare package Medium = MediumWater) annotation (Placement(transformation(extent={{14,-42},{34,-22}})));
   
  .FuelCell.Sensors.WaterMultiDisplaySensor multiDisplaySensor4(redeclare package Medium = MediumWater) annotation (Placement(transformation(extent={{-36,-50},{-16,-30}})));
  .FuelCell.Sensors.WaterMultiDisplaySensor multiDisplaySensor1(redeclare package Medium = MediumWater) annotation (Placement(transformation(extent={{14,-50},{34,-30}})));
  .FuelCell.Sensors.GasMultiDisplaySensor gasSensor3(redeclare package Medium = Medium_an) annotation (Placement(transformation(extent={{-78,20},{-58,40}})));
  .FuelCell.Sensors.GasMultiDisplaySensor gasSensor2(redeclare package Medium = Medium_an) annotation (Placement(transformation(extent={{24,20},{44,40}})));
  .FuelCell.Sensors.GasMultiDisplaySensor gasSensor4(redeclare package Medium = Medium_cath) annotation (Placement(transformation(extent={{-78,-30},{-58,-10}})));
  .FuelCell.Sensors.GasMultiDisplaySensor gasSensor1(redeclare package Medium = Medium_cath) annotation (Placement(transformation(extent={{26,-30},{46,-10}})));
  .Modelica.Blocks.Sources.Ramp ramp(
    startTime = 3600,
    offset = m_flow_water_min,
    duration = 0,
    height = m_flow_water_max - m_flow_water_min) annotation(Placement(transformation(extent = {{-81.91734927990456,-47.917349279904556},{-70.08265072009544,-36.082650720095444}},origin = {0.0,0.0},rotation = 0.0)));
  .Modelica.Blocks.Sources.Ramp ramp2(
    height = if constant_current then 0 else m_flow_an_max - m_flow_an_min, 
    duration = if constant_current then 0 else simulation_time,
    offset = m_flow_an_min,
    startTime = if constant_current then simulation_time else 0) annotation(Placement(transformation(extent = {{-125.91734927990456,40.082650720095444},{-114.08265072009544,51.917349279904556}},origin = {0.0,0.0},rotation = 0.0)));
  .Modelica.Blocks.Sources.Ramp ramp3(
    height = if constant_current then 0 else m_flow_cath_max - m_flow_cath_min,
    duration = if constant_current then 0 else simulation_time,
    offset = m_flow_cath_min,
    startTime = if constant_current then simulation_time else 0) annotation(Placement(transformation(extent = {{-126.18678582346512,-9.776568124258443},{-114.352087263656,2.0581304355506695}},origin = {0.0,0.0},rotation = 0.0)));
   
          
equation       
        
  connect(current.n, ground. p) annotation (Line(
      points={{7.000000000000001,64},{7.000000000000001,71.87114738360049},{25.040856259459844,71.87114738360049}},
      color={0,0,255},
      smooth=Smooth.None));
  connect(current.p,coolStack.pin_p) annotation(Line(points = {{-7.000000000000001,64},{-10,64},{-10,36.5}},color = {0,0,255}));
  connect(current.n,coolStack.pin_n) annotation(Line(points = {{7.000000000000001,64},{16.362207960984374,64},{16.362207960984374,36.5},{6,36.5}},color = {0,0,255}));

  connect(multiDisplaySensor4.u, display_phTmdot_w_in. y) annotation (Line(points={{-25.95,-39.95},{-25.547635814819174,-32}}, color={0,0,0}));
  connect(multiDisplaySensor1.u, display_phTmdot_w_out. y) annotation (Line(points={{24.05,
          -39.95},{24.05,-35.975},{24,-35.975},{24,-32}}, color={0,0,0}));
  connect(multiDisplaySensor4.portA, sourceW.fluidPort) annotation (Line(points={{-32,-40},{-39,-40}}, color={0,0,255}));
  connect(multiDisplaySensor4.portB, coolStack.feed_cooling) annotation (Line(points={{-20,-40},{-10,-40},{-10,3.5}}, color={0,0,255}));
  connect(multiDisplaySensor1.portB, sinkP.fluidPort) annotation (Line(points={{30,-40},{39,-40}}, color={0,0,255}));
  connect(multiDisplaySensor1.portA, coolStack.drain_cooling) annotation (Line(points={{18,-40},{6,-40},{6,3.5}}, color={0,0,255}));
  connect(display_phTmdot_an_in.y, gasSensor3. u) annotation (Line(points={{-68,40},{
          -68,30.05},{-67.95,30.05}},
                                  color={0,0,0}));
  connect(display_an_in.data, gasSensor3. u) annotation (Line(points={{-40.462753249644805,37.54955698575693},{-56,37.54955698575693},{-56,30.05},{-67.95,30.05}},
                                                         color={0,0,0}));
  connect(gasSensor3.portA, flowAnode.fluidPort) annotation (Line(points={{-74,30},{-82.93906632502666,30},{-82.93906632502666,29.768425703538625}}, color={255,128,0}));
  connect(gasSensor3.portB, coolStack.feed_an) annotation (Line(points={{-62,30},{-20,30}}, color={255,128,0}));
  connect(gasSensor2.u, display_phTmdot_an_out. y) annotation (Line(points={{34.05,
          30.05},{34.05,35.025},{34,35.025},{34,40}},
                                               color={0,0,0}));
  connect(gasSensor2.u, display_an_out.data) annotation (Line(points={{34.05,30.05},{34.05,32},{62,32},{62,37.2941}},
                                                   color={0,0,0}));
  connect(gasSensor2.portA, coolStack.drain_an) annotation (Line(points={{28,30},{16,30}}, color={255,128,0}));
  connect(gasSensor2.portB, sinkAnode.fluidPort) annotation (Line(points={{40,30},{80.22993225215629,30},{80.22993225215629,29.6149661260782}}, color={255,128,0}));
  connect(gasSensor4.u, display_phTmdot_cath_in. y) annotation (Line(points={{-67.95,-19.95},{-67.95,-14.975},{-68,-14.975},{-68,-10}},
                                                     color={0,0,0}));
  connect(gasSensor4.u, display_cath_in. data) annotation (Line(points={{-67.95,-19.95},{-42.744560661301904,-19.95},{-42.744560661301904,-14}}, color={0,0,0}));
  connect(gasSensor4.portA, flowCathode.fluidPort) annotation (Line(points={{-74,-20},{-83,-20}}, color={255,128,0}));
  connect(gasSensor4.portB, coolStack.feed_cath) annotation (Line(points={{-62,-20},{-20,-20},{-20,10}}, color={255,128,0}));
  connect(gasSensor1.u, display_cath_out.data) annotation (Line(points={{36.05,-19.95},{62.00000000000004,-19.95},{62.00000000000004,-13}},    color={0,0,0}));
  connect(sinkCathode.fluidPort, gasSensor1.portB) annotation (Line(points={{81,-20},{42,-20}}, color={255,128,0}));
  connect(gasSensor1.portA, coolStack.drain_cath) annotation (Line(points={{30,-20},{16,-20},{16,10}}, color={255,128,0}));
  connect(gasSensor1.u, display_phTmdot_cath_out.y) annotation (Line(points={{36.05,
          -19.95},{36.05,-14.975},{36,-14.975},{36,-10}}, color={0,0,0}));
  connect(ramp.y,sourceW.m_flow_in) annotation(Line(points = {{-69.49091579210499,-42},{-61.95890105565714,-42},{-61.95890105565714,-30},{-54,-30}},color = {0,0,127}));
  connect(ramp2.y,flowAnode.m_flow_in) annotation(Line(points = {{-113.49091579210499,46},{-97.93906632502666,46},{-97.93906632502666,39.768425703538625}},color = {0,0,127}));
  connect(ramp3.y,flowCathode.m_flow_in) annotation(Line(points = {{-113.76035233566554,-3.8592188443538866},{-98,-3.8592188443538866},{-98,-10}},color = {0,0,127}));
    
  annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-120,
            -100},{120,100}}), graphics={
        Rectangle(
          extent={{-96,-60},{0,-98}},
          lineColor={215,215,215},
          fillColor={215,215,215},
          fillPattern=FillPattern.Solid,
          radius=2),
        Text(
          extent={{-88,-66},{-70,-72}},
          lineColor={0,0,0},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid,
          textString="Voltage [V]"),
        Text(
          extent={{-66,-66},{-48,-72}},
          lineColor={0,0,0},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid,
          textString="Power [W]"),
        Line(points={{-96,-66},{0,-66}},  color={0,0,0}),
        Text(
          extent={{-94,-66},{-68,-60}},
          lineColor={0,0,0},
          textStyle={TextStyle.Bold},
          textString="Stack"),
        Text(
          extent={{-88,-82},{-70,-88}},
          lineColor={0,0,0},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid,
          textString="Current [A]"),
        Text(
          extent={{-46,-64},{-6,-74}},
          lineColor={0,0,0},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid,
          textString="Current density [mA/cm^2]"),
        Text(
          extent={{-38,-82},{-16,-88}},
          lineColor={0,0,0},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid,
          textString="Cooling [W]"),
        Text(
          extent={{-66,-82},{-44,-88}},
          lineColor={0,0,0},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid,
          textString="T cell [C]")}),
    experiment(StopTime=1000),
    __Dymola_experimentSetupOutput(equdistant=true),
    Documentation(revisions="<html>
Copyright &copy; 2004-2024, MODELON AB <br /> The use of this software component is regulated by the licensing conditions for Modelon Libraries. <br />This copyright notice must, unaltered, accompany all components that are derived from, copied from, <br />or by other means have their origin from any Modelon Library.
</html>", info="<html>
<h4>TestPEMFCCoolStack</h4>
<p>Experiment with a proton exchange membrane fuel cell with cooling. All boundary conditions are constant, making the model converge to steady state heat production and temperature.</p>
</html>"),
    Icon(coordinateSystem(extent={{-120,-100},{120,100}})));

end CoolStack_oc;
