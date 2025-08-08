within PEMFCModel.Stacks.Experiments.PEMFC;

model CoolStack_constantpressure_humidified "Cooled PEMFC stack test model"
    
  extends .Modelon.Icons.Experiment;

  replaceable package Medium_an = .PEMFCModel4par.Medium.NASASteamH2N2O2 constrainedby
    .FuelCell.Media.Templates.ReactionGas "Anodic medium";
  replaceable package Medium_cath = .FuelCell.Media.PreDefined.IdealGases.NASAMoistAir constrainedby
    .FuelCell.Media.Templates.ReactionGas "Cathodic medium";
  package MediumWater = .FuelCell.Media.PreDefined.TwoPhase.WaterIF97 "Water medium";
  
  parameter Real V_flow_an_min(unit = "NLPM") = 29.59 "Minimum value of anodic dry gas inlet volumetric flow rate at variable load" annotation(Dialog(group = "Anode side", enable = not constant_current));
  parameter Real V_flow_an_max(unit = "NLPM") = 73.974 "Maximum value of anodic dry gas inlet volumetric flow rate at variable load" annotation(Dialog(group = "Anode side", enable = not constant_current));
  parameter Real V_flow_an(unit = "NLPM") = 73.974 "Anodic dry gas inlet volumetric flow rate at constant load" annotation(Dialog(group = "Anode side", enable = constant_current));
  parameter .Modelica.Units.SI.Pressure p_an_in = 2.2e5 "Start inlet pressure at anode side" annotation(Dialog(group = "Anode side"));
  parameter .Modelica.Units.SI.Pressure p_an_out = p_an_in - 0.1e5 "Outlet pressure at anode side" annotation(Dialog(group = "Anode side"));
  parameter .Modelica.Units.SI.Temperature T_an_in = 70+273.15 "Anodic gas inlet temperature" annotation(Dialog(group = "Anode side"));
  parameter .Modelica.Units.SI.MoleFraction y_dry_an[Medium_an.nS] = {0.7,0,0.7,0} "Anodic dry gas inlet molar-based composition"  annotation(Dialog(group = "Anode side"));
  parameter .Modelica.Units.SI.Temperature T_dew_an = 52.5+273.15 "Anodic gas dew point" annotation(Dialog(group = "Anode side"));
        
  parameter Real V_flow_cath_min(unit = "NLPM") = 59.23 "Minimum value of cathodic dry gas inlet volumetric flow rate at variable load" annotation(Dialog(group = "Cathode side", enable = not constant_current));
  parameter Real V_flow_cath_max(unit = "NLPM") = 148.076 "Maximum value of cathodic dry gas inlet volumetric flow rate at variable load" annotation(Dialog(group = "Cathode side", enable = not constant_current));
  parameter Real V_flow_cath(unit = "NLPM") = 148.076 "Cathodic dry gas inlet volumetric flow rate at constant load" annotation(Dialog(group = "Cathode side", enable = constant_current));
  parameter .Modelica.Units.SI.Pressure p_cath_in = 2e5 "Start inlet pressure at cathode side" annotation(Dialog(group = "Cathode side"));
  parameter .Modelica.Units.SI.Pressure p_cath_out = p_cath_in - 0.1e5 "Outlet pressure at cathode side" annotation(Dialog(group = "Cathode side"));
  parameter .Modelica.Units.SI.Temperature T_cath_in = 70+273.15 "Cathodic gas inlet temperature" annotation(Dialog(group = "Cathode side"));
  parameter .Modelica.Units.SI.MassFraction y_dry_cath[Medium_cath.nS] = {0,0,0,0.79,0.21} "Cathodic dry gas inlet molar-based composition" annotation(Dialog(group = "Cathode side"));
  parameter .Modelica.Units.SI.Temperature T_dew_cath = 52.5+273.15 "Cathodic gas dew point" annotation(Dialog(group = "Cathode side"));

  parameter Real V_flow_water_min(unit = "NLPM") = 4.95 "Minimum value of water inlet volumetric flow rate at variable load" annotation(Dialog(group = " Cooling", enable = not constant_current));
  parameter Real V_flow_water_max(unit = "NLPM") = 4.95 "Maximum value of water inlet volumetric flow rate at variable load" annotation(Dialog(group = " Cooling", enable = not constant_current));
  parameter Real V_flow_water(unit = "NLPM") = 4.95 "Water inlet volumetric flow rate at constant load" annotation(Dialog(group = " Cooling", enable = constant_current));      
  constant .Modelica.Units.SI.Density rho_water = 1000 "Constant water density";
  parameter .Modelica.Units.SI.MassFlowRate m_flow_water_min = V_flow_water_min/1000/60*rho_water  annotation(Dialog(group = " Cooling", enable = false));
  parameter .Modelica.Units.SI.MassFlowRate m_flow_water_max = V_flow_water_max/1000/60*rho_water  annotation(Dialog(group = " Cooling", enable = false));
  parameter .Modelica.Units.SI.Temperature T_water_in = 65+273.15 "Water inlet temperature" annotation(Dialog(group = " Cooling"));
  parameter .Modelica.Units.SI.Pressure p_water_in = 2.4e5 "Start water inlet pressure" annotation(Dialog(group = " Cooling"));
  parameter .Modelica.Units.SI.Pressure p_water_out = p_water_in - 0.4e5 "Water outlet pressure" annotation(Dialog(group = " Cooling"));
  //parameter .Modelica.Units.SI.MassFlowRate m_out=m_air*X_air[cath[1]]*(1 - 1/cath_stoich)
    //"calculated mass flow rate of oxygen gas at cathode outlet, may be used for control of the inlet cathode flow rate";
  parameter .Modelica.Units.SI.Current I_min = 180 "Minimum value of current at variable load" annotation(Dialog(group = "Load", enable = not constant_current)); 
  parameter .Modelica.Units.SI.Current I_max = 450 "Maximum value of current at variable load" annotation(Dialog(group = "Load", enable = not constant_current));
  parameter .Modelica.Units.SI.Current I = 450 "Current at constant load" annotation(Dialog(group = "Load", enable = constant_current));
  parameter Boolean constant_current = true "If true, load is constant" annotation(Dialog(group = "Load"));
 
  parameter .Modelica.Units.SI.Time simulation_time = 3600 "Simulation time" annotation(Dialog(tab = "Simulation options"));         
  
  .PEMFCModel.Humidification humidification_an(
    redeclare package Medium = Medium_an,
    T_in = T_an_in ,
    y_dry_in = y_dry_an,
    T_dew = T_dew_an) annotation(Placement(transformation(extent = {{-130.3122469373013,23.687753062698697},{-117.6877530626987,36.3122469373013}},origin = {0.0,0.0},rotation = 0.0)));
  .PEMFCModel.Humidification humidification_cath(
    redeclare package Medium = Medium_cath,
    T_in = T_cath_in ,
    y_dry_in = y_dry_cath,
    T_dew = T_dew_cath) annotation(Placement(transformation(extent = {{-130.3122469373013,-24.312246937301303},{-117.6877530626987,-11.687753062698697}},origin = {0.0,0.0},rotation = 0.0)));
           
  .FuelCell.Sources.GasFlowBoundary flowAnode(
    redeclare package Medium = Medium_an,
    T = T_an_in,
    use_flow_in = true,
    use_Th_in = false,use_X_in = true) annotation (Placement(transformation(
          extent={{-102.19817221640716,19.768425703538625},{-82.19817221640716,39.768425703538625}},
                                     rotation=0.0,origin = {0.0,0.0})));

  .FuelCell.Sources.GasPressureBoundary sinkAnode(
    
    T = T_an_in,
    redeclare package Medium = Medium_an,use_p_in = false,use_X_in = true,p = p_an_out,X = Medium_an.moleToMassFractions(y_dry_an,Medium_an.MMX),use_Th_in = false) annotation (Placement(transformation(extent={{99.22993225215629,19.6149661260782},{79.22993225215629,39.6149661260782}},
          rotation=0.0,origin = {0.0,0.0})));

  .Modelica.Electrical.Analog.Basic.Ground ground annotation (Placement(transformation(extent={{19.040856259459844,59.87114738360049},{31.040856259459844,71.87114738360049}},   rotation=
           0.0,origin = {0.0,0.0})));
        
  .Modelica.Electrical.Analog.Sources.RampCurrent current(
    startTime = if constant_current then simulation_time else 0,
    duration = if constant_current then 0 else simulation_time,
    offset = if constant_current then I else I_min,
    I = if constant_current then I else I_max - I_min) annotation (Placement(transformation(
        origin={0.1640169024502267,64.78525109495013},
        extent={{7.0,-7.0},{-7.0,7.0}},
        rotation=-180.0)));
        
  .FuelCell.Sources.GasFlowBoundary flowCathode(
    redeclare package Medium = Medium_cath,
    T = T_cath_in,
    use_flow_in = true,
    use_Th_in = false,use_X_in = true) annotation (Placement(transformation(
          extent={{-102.0,-30.0},{-82.0,-10.0}},
                                      rotation=0.0,origin = {0.0,0.0})));
        
  .FuelCell.Sources.GasPressureBoundary sinkCathode(
    redeclare package Medium = Medium_cath,
    T = T_cath_in,use_p_in = false,use_Th_in = false,use_X_in = true,p = p_cath_out,X = Medium_cath.moleToMassFractions(y_dry_cath,Medium_cath.MMX)) annotation (Placement(transformation(extent={{99.61496612607814,-30.38503387392182},{79.61496612607814,-10.38503387392182}},
          rotation=0.0,origin = {0.0,0.0})));
        
  .PEMFCModel.Stacks.PEMFC.Empirical.CoolStack coolStack(
    redeclare package Medium_an = Medium_an,
    redeclare package Medium_cath = Medium_cath,
    redeclare package Medium_cooling = MediumWater,
    X_feed_an = humidification_an.x_wet_in,
    X_feed_cath = humidification_cath.x_wet_in,p_start_out_anode = sinkAnode.p,p_start_in_anode = p_an_in,m_flow_start_cathode = 1e-3,p_start_out_cathode = sinkCathode.p,p_start_in_cathode = p_cath_in,T_start_in_anode = flowAnode.T,T_start_in_cathode = flowCathode.T,T_start_out_cathode = flowCathode.T,T_start_out_anode = flowAnode.T,positiveFlow_cooling = true,from_dp_cooling = false,positiveFlow_anode = true,positiveFlow_cathode = true,from_dp_anode = false,from_dp_cathode = false,p_start_in_cooling = p_water_in,T_start_in_cooling = sourceW.T,T_start_out_cooling = sourceW.T,m_flow_start_cooling = m_flow_water_min, h_inflow_an = flowAnode.fluidPort.h_outflow, h_inflow_cath = flowCathode.fluidPort.h_outflow, h_inflow_cooling = sourceW.fluidPort.h_outflow,p_start_out_cooling = sinkP.p,enable_setting = true,dp_smooth = 1e-2,mflow_smooth = 1e-8,dp_smooth_cooling = 1e-2,mflow_smooth_cooling = 1e-8,redeclare replaceable model Friction_anode = .Modelon.ThermoFluid.FlowChannels.PipeResistances.SinglePhase.LinearOperatingPointLoss(d0 = 0.4,dp0 = 0.1e5,m_flow0 = 1e-4),redeclare replaceable model Friction_cathode = .Modelon.ThermoFluid.FlowChannels.PipeResistances.SinglePhase.LinearOperatingPointLoss(d0 = 1.2,dp0 = 0.15e5,m_flow0 = 1e-3),redeclare replaceable model HeatTransfer_anode = .Modelon.ThermoFluid.FlowChannels.HeatTransfer.SinglePhase.SinglePhase,redeclare replaceable model HeatTransfer_cathode = .Modelon.ThermoFluid.FlowChannels.HeatTransfer.SinglePhase.SinglePhase,redeclare replaceable model Friction_cooling = .Modelon.ThermoFluid.FlowChannels.PipeResistances.TwoPhase.DensityProfileFriction(mflow0 = m_flow_water_min,p0_in = p_water_in,p0_out = p_water_out),redeclare replaceable model HeatTransfer_cooling = .Modelon.ThermoFluid.FlowChannels.HeatTransfer.TwoPhase.SimpleTwoPhase,N = 5,A_cell = 300e-4,L_anode = 0.0813,D_anode = 0.023,D_cathode = 0.0445,L_cathode = 0.0813,L_cooling = 0.0813,D_cooling = 0.0375,diffusiveSpecies = {"O2","N2"},n_cell = 11,M_stack = 10.1,z = 75e-6,j_0 = 10,m_conc = 3e-4,alpha = 0.5,c1 = 5.75e-3,c_stack = 10e8,m_flow_start_anode = 1e-4,n_conc = 3e-4) annotation (Placement(transformation(extent={{-22.0,-2.0},{18.0,38.0}},rotation = 0.0,origin = {0.0,0.0})));

   .FuelCell.Sources.WaterFlowBoundary sourceW(redeclare package Medium = MediumWater, T = T_water_in,use_flow_in = true,use_Th_in = false) annotation (Placement(transformation(extent={{-57.62305377241874,-50.0},{-37.62305377241874,-30.0}},rotation = 0.0,origin = {0.0,0.0})));
        
  .FuelCell.Sources.WaterPressureBoundary sinkP(redeclare package Medium = MediumWater, T = T_water_in,use_p_in = false,p = p_water_out) annotation (Placement(transformation(extent={{-10.0,10.0},{10.0,-10.0}},
        rotation=180.0,
        origin={48.58183823418737,-40.94544117720845})));
        
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
  
  .Modelica.Blocks.Sources.Ramp ramp_m_flow_water(
    startTime = 3600,
    offset = m_flow_water_min,
    duration = 0,
    height = m_flow_water_max - m_flow_water_min) annotation(Placement(transformation(extent = {{-82.70681533891384,-48.461718351330504},{-70.87211677910473,-36.62701979152139}},origin = {0.0,0.0},rotation = 0.0)));
  .Modelica.Blocks.Sources.Ramp ramp_V_flow_cath_dry(
    height = if constant_current then 0 else V_flow_cath_max - V_flow_cath_min,
    duration = if constant_current then 0 else simulation_time,
    offset = if constant_current then V_flow_cath else V_flow_cath_min,
    startTime = if constant_current then simulation_time else 0) annotation(Placement(transformation(extent = {{-167.91734927990456,-37.917349279904556},{-156.08265072009544,-26.082650720095444}},origin = {0.0,0.0},rotation = 0.0)));
  .Modelica.Blocks.Sources.Ramp ramp_V_flow_an_dry(
    height = if constant_current then 0 else V_flow_an_max - V_flow_an_min,
    duration = if constant_current then 0 else simulation_time,
    offset = if constant_current then V_flow_an else V_flow_an_min,
    startTime = if constant_current then simulation_time else 0) annotation(Placement(transformation(extent = {{-169.91734927990456,14.082650720095444},{-158.08265072009544,25.917349279904556}},origin = {0.0,0.0},rotation = 0.0)));
    .Modelica.Blocks.Sources.RealExpression p_an_inlet(y = p_an_in) annotation(Placement(transformation(extent = {{-174.30768824601444,31.421598899305913},{-154.30768824601444,51.42159889930591}},origin = {0.0,0.0},rotation = 0.0)));
    .Modelica.Blocks.Sources.RealExpression p_cath_inlet(y = p_cath_in) annotation(Placement(transformation(extent = {{-172.0,-18.0},{-152.0,2.0}},origin = {0.0,0.0},rotation = 0.0)));
           
equation       
        
  connect(current.n, ground. p) annotation (Line(
      points={{7.164016902450228,64.78525109495013},{7.164016902450228,71.87114738360049},{25.040856259459844,71.87114738360049}},
      color={0,0,255},
      smooth=Smooth.None));
  connect(current.p,coolStack.pin_p) annotation(Line(points = {{-6.835983097549774,64.78525109495013},{-10,64.78525109495013},{-10,34.5}},color = {0,0,255}));
  connect(current.n,coolStack.pin_n) annotation(Line(points = {{7.164016902450228,64.78525109495013},{16.362207960984374,64.78525109495013},{16.362207960984374,34.5},{6,34.5}},color = {0,0,255}));

  connect(multiDisplaySensor4.u, display_phTmdot_w_in. y) annotation (Line(points={{-25.95,-39.95},{-25.547635814819174,-32}}, color={0,0,0}));
  connect(multiDisplaySensor1.u, display_phTmdot_w_out. y) annotation (Line(points={{24.05,
          -39.95},{24.05,-35.975},{24,-35.975},{24,-32}}, color={0,0,0}));
  connect(multiDisplaySensor4.portA, sourceW.fluidPort) annotation (Line(points={{-32,-40},{-38.62305377241874,-40}}, color={0,0,255}));
  connect(multiDisplaySensor4.portB, coolStack.feed_cooling) annotation (Line(points={{-20,-40},{-10,-40},{-10,1.5}}, color={0,0,255}));
  connect(multiDisplaySensor1.portB, sinkP.fluidPort) annotation (Line(points={{30,-40},{39.58183823418737,-40},{39.58183823418737,-40.94544117720845}}, color={0,0,255}));
  connect(multiDisplaySensor1.portA, coolStack.drain_cooling) annotation (Line(points={{18,-40},{6,-40},{6,1.5}}, color={0,0,255}));
  connect(display_phTmdot_an_in.y, gasSensor3. u) annotation (Line(points={{-68,40},{
          -68,30.05},{-67.95,30.05}},
                                  color={0,0,0}));
  connect(display_an_in.data, gasSensor3. u) annotation (Line(points={{-40.462753249644805,37.54955698575693},{-56,37.54955698575693},{-56,30.05},{-67.95,30.05}},
                                                         color={0,0,0}));
  connect(gasSensor3.portA, flowAnode.fluidPort) annotation (Line(points={{-74,30},{-83.19817221640716,30},{-83.19817221640716,29.768425703538625}}, color={255,128,0}));
  connect(gasSensor3.portB, coolStack.feed_an) annotation (Line(points={{-62,30},{-20,30},{-20,28}}, color={255,128,0}));
  connect(gasSensor2.u, display_phTmdot_an_out. y) annotation (Line(points={{34.05,
          30.05},{34.05,35.025},{34,35.025},{34,40}},
                                               color={0,0,0}));
  connect(gasSensor2.u, display_an_out.data) annotation (Line(points={{34.05,30.05},{34.05,32},{62,32},{62,37.2941}},
                                                   color={0,0,0}));
  connect(gasSensor2.portA, coolStack.drain_an) annotation (Line(points={{28,30},{16,30},{16,28}}, color={255,128,0}));
  connect(gasSensor2.portB, sinkAnode.fluidPort) annotation (Line(points={{40,30},{80.22993225215629,30},{80.22993225215629,29.6149661260782}}, color={255,128,0}));
  connect(gasSensor4.u, display_phTmdot_cath_in. y) annotation (Line(points={{-67.95,-19.95},{-67.95,-14.975},{-68,-14.975},{-68,-10}},
                                                     color={0,0,0}));
  connect(gasSensor4.u, display_cath_in. data) annotation (Line(points={{-67.95,-19.95},{-42.744560661301904,-19.95},{-42.744560661301904,-14}}, color={0,0,0}));
  connect(gasSensor4.portA, flowCathode.fluidPort) annotation (Line(points={{-74,-20},{-83,-20}}, color={255,128,0}));
  connect(gasSensor4.portB, coolStack.feed_cath) annotation (Line(points={{-62,-20},{-20,-20},{-20,8}}, color={255,128,0}));
  connect(gasSensor1.u, display_cath_out.data) annotation (Line(points={{36.05,-19.95},{62.00000000000004,-19.95},{62.00000000000004,-13}},    color={0,0,0}));
  connect(sinkCathode.fluidPort, gasSensor1.portB) annotation (Line(points={{80.61496612607814,-20.38503387392182},{80.61496612607814,-20},{42,-20}}, color={255,128,0}));
  connect(gasSensor1.portA, coolStack.drain_cath) annotation (Line(points={{30,-20},{16,-20},{16,8}}, color={255,128,0}));
  connect(gasSensor1.u, display_phTmdot_cath_out.y) annotation (Line(points={{36.05,
          -19.95},{36.05,-14.975},{36,-14.975},{36,-10}}, color={0,0,0}));
  connect(humidification_an.m_flow_wet_in,flowAnode.m_flow_in) annotation(Line(points = {{-116.40952305789519,32.588021244293536},{-110.40952305789519,32.588021244293536},{-110.40952305789519,45.768425703538625},{-98.19817221640716,45.768425703538625},{-98.19817221640716,39.768425703538625}},color = {0,0,127}));
  connect(humidification_an.x_wet_in,flowAnode.X_in) annotation(Line(points = {{-116.18859441508964,27.91695851069057},{-110,27.91695851069057},{-110,45.768425703538625},{-86.19817221640716,45.768425703538625},{-86.19817221640716,39.768425703538625}},color = {0,0,127}));
  connect(humidification_cath.m_flow_wet_in,flowCathode.m_flow_in) annotation(Line(points = {{-116.40952305789519,-15.411978755706466},{-110.40952305789519,-15.411978755706466},{-110.40952305789519,-4},{-98,-4},{-98,-10}},color = {0,0,127}));
  connect(humidification_cath.x_wet_in,flowCathode.X_in) annotation(Line(points = {{-116.18859441508964,-20.08304148930943},{-110,-20.08304148930943},{-110,-4},{-86,-4},{-86,-10}},color = {0,0,127}));
  connect(ramp_m_flow_water.y,sourceW.m_flow_in) annotation(Line(points = {{-70.28038185111427,-42.54436907142595},{-61.95890105565714,-42.54436907142595},{-61.95890105565714,-30},{-53.62305377241874,-30}},color = {0,0,127}));
  connect(humidification_an.x_wet_in,sinkAnode.X_in) annotation(Line(points = {{-116.18859441508964,27.91695851069057},{-110.18859441508964,27.91695851069057},{-110.18859441508964,84},{83.22993225215629,84},{83.22993225215629,39.6149661260782}},color = {0,0,127}));
  connect(humidification_cath.x_wet_in,sinkCathode.X_in) annotation(Line(points = {{-116.18859441508964,-20.08304148930943},{-110,-20.08304148930943},{-110,18},{83.61496612607814,18},{83.61496612607814,-10.38503387392182}},color = {0,0,127}));
    connect(p_cath_inlet.y,humidification_cath.p_in) annotation(Line(points = {{-151,-8},{-141.2873481623808,-8},{-141.2873481623808,-14.212651837619218},{-131.57469632476156,-14.212651837619218}},color = {0,0,127}));
    connect(p_an_inlet.y,humidification_an.p_in) annotation(Line(points = {{-153.30768824601444,41.42159889930591},{-142.2873481623808,41.42159889930591},{-142.2873481623808,33.78734816238078},{-131.57469632476156,33.78734816238078}},color = {0,0,127}));
    connect(ramp_V_flow_an_dry.y,humidification_an.V_flow_dry_in) annotation(Line(points = {{-157.490915792105,20},{-144.53280605843327,20},{-144.53280605843327,26.21265183761922},{-131.57469632476156,26.21265183761922}},color = {0,0,127}));
    connect(ramp_V_flow_cath_dry.y,humidification_cath.V_flow_dry_in) annotation(Line(points = {{-155.490915792105,-32},{-143.53280605843327,-32},{-143.53280605843327,-21.78734816238078},{-131.57469632476156,-21.78734816238078}},color = {0,0,127}));
    
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

end CoolStack_constantpressure_humidified;
