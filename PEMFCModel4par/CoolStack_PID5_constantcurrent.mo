within PEMFCModel4par;

model CoolStack_PID5_constantcurrent "Cooled PEMFC stack test model"
    
  extends .Modelon.Icons.Experiment;

  replaceable package Medium_an = .PEMFCModel4par.Medium.NASASteamH2N2O2 constrainedby
    .FuelCell.Media.Templates.ReactionGas "Anodic medium";
  replaceable package Medium_cath = .FuelCell.Media.PreDefined.IdealGases.NASAMoistAir constrainedby
    .FuelCell.Media.Templates.ReactionGas "Cathodic medium";
  package MediumWater = .FuelCell.Media.PreDefined.TwoPhase.WaterIF97 "Water medium";
   
  parameter Real V_flow_an(unit = "NLPM") = 29.59 "Anodic dry gas inlet volumetric flow rate" annotation(Dialog(group = "Anode side"));
  parameter .Modelica.Units.SI.Pressure p_an_in = 1.9e5 "Anodic gas inlet pressure" annotation(Dialog(group = "Anode side"));      
  parameter .Modelica.Units.SI.Pressure dp_an = 0.04e5 "Pressure drop at anode side" annotation(Dialog(group = "Anode side"));
  parameter .Modelica.Units.SI.Temperature T_an_in = 70+273.15 "Anodic gas inlet temperature" annotation(Dialog(group = "Anode side"));
  parameter .Modelica.Units.SI.MoleFraction y_dry_an[Medium_an.nS] = {0.7,0,0.7,0} "Anodic dry gas inlet molar-based composition"  annotation(Dialog(group = "Anode side"));
  parameter .Modelica.Units.SI.Temperature T_dew_an = 52.5+273.15 "Anodic gas dew point" annotation(Dialog(group = "Anode side"));
        
  parameter Real V_flow_cath(unit = "NLPM") = 59.23 "Cathodic dry gas inlet volumetric flow rate" annotation(Dialog(group = "Cathode side"));
  parameter .Modelica.Units.SI.Pressure p_cath_in = 1.7e5 "Cathodic gas inlet pressure" annotation(Dialog(group = "Cathode side"));  
  parameter Real dp_cath = 0.1e5 "Pressure drop at cathode side" annotation(Dialog(group = "Cathode side"));
  parameter .Modelica.Units.SI.Temperature T_cath_in = 70+273.15 "Cathodic gas inlet temperature" annotation(Dialog(group = "Cathode side"));
  parameter .Modelica.Units.SI.MassFraction y_dry_cath[Medium_cath.nS] = {0,0,0,0.79,0.21} "Cathodic dry gas inlet molar-based composition" annotation(Dialog(group = "Cathode side"));
  parameter .Modelica.Units.SI.Temperature T_dew_cath = 52.5+273.15 "Cathodic gas dew point" annotation(Dialog(group = "Cathode side"));

  parameter Real V_flow_water(unit = "LPM") = 4.95 "Water inlet volumetric flow rate" annotation(Dialog(group = " Cooling"));      
  constant .Modelica.Units.SI.Density rho_water = 1000 "Water density";
  parameter .Modelica.Units.SI.MassFlowRate m_flow_water = V_flow_water/1000/60*rho_water  annotation(Dialog(group = " Cooling", enable = false));
  parameter .Modelica.Units.SI.Pressure p_water_in = 2.1e5 "Water inlet pressure" annotation(Dialog(group = " Cooling"));
  parameter .Modelica.Units.SI.Temperature T_water_in = 65+273.15 "Water inlet temperature" annotation(Dialog(group = " Cooling"));
  
  //parameter .Modelica.Units.SI.MassFlowRate m_out=m_air*X_air[cath[1]]*(1 - 1/cath_stoich)
    //"calculated mass flow rate of oxygen gas at cathode outlet, may be used for control of the inlet cathode flow rate";
  parameter .Modelica.Units.SI.Current I = 180 "Load current" annotation(Dialog(group = "Load"));
  parameter .Modelica.Units.SI.Voltage V_exp = 0.742 "Stack voltage from experiments" annotation(Dialog(group = "Load"));
 
  parameter Real m_flow_an_nom = 6.072e-4 "Nominal value of anodic wet gas inlet mass flow rate for linear friction loss model" annotation(Dialog(tab = "Nominal conditions", group = "Anode side"));
  parameter Real m_flow_cath_nom = 3.327e-3 "Nominal value of cathodic wet gas inlet mass flow rate for linear friction loss model" annotation(Dialog(tab = "Nominal conditions", group = "Cathode side"));
  parameter Real m_flow_water_nom = 8.25e-2 "Nominal value of water inlet mass flow rate for density profile friction model" annotation(Dialog(tab = "Nominal conditions", group = "Cooling"));
  parameter Real rho_an_nom = 0.80 "Nominal value of anodic wet gas inlet density for linear friction loss model"  annotation(Dialog(tab = "Nominal conditions", group = "Anode side"));
  parameter Real rho_cath_nom = 1.98 "Nominal value of cathodic wet gas inlet density for linear friction loss model"  annotation(Dialog(tab = "Nominal conditions", group = "Cathode side"));
  parameter Real dp_an_nom = 0.20e5 "Nominal value of pressure drop at anode for linear friction loss model" annotation(Dialog(tab = "Nominal conditions", group = "Anode side"));
  parameter Real dp_cath_nom = 0.225e5 "Nominal value of pressure drop at cathode for linear friction loss model" annotation(Dialog(tab = "Nominal conditions", group = "Cathode side"));
  parameter Real p_water_in_nom = 2.4e5 "Nominal value of water inlet pressure for density profile friction model" annotation(Dialog(tab = "Nominal conditions", group = "Cooling"));
  parameter Real dp_water_nom = 0.43e5 "Nominal value of cooling channel pressure drop for density profile friction model" annotation(Dialog(tab = "Nominal conditions", group = "Cooling"));
        
  Real err_rel(unit = "%") = abs(coolStack.subStack.cell.V_cell - V_exp)/V_exp*100;
  
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
          extent={{-101.93906632502666,19.768425703538625},{-81.93906632502666,39.768425703538625}},
                                     rotation=0.0,origin = {0.0,0.0})));

  .FuelCell.Sources.GasPressureBoundary sinkAnode(
    
    T = T_an_in,
    redeclare package Medium = Medium_an,use_p_in = true,use_X_in = true) annotation (Placement(transformation(extent={{99.22993225215629,19.6149661260782},{79.22993225215629,39.6149661260782}},
          rotation=0.0,origin = {0.0,0.0})));

  .Modelica.Electrical.Analog.Basic.Ground ground annotation (Placement(transformation(extent={{24.0,52.0},{36.0,64.0}},   rotation=
           0.0,origin = {0.0,0.0})));
        
  .FuelCell.Sources.GasFlowBoundary flowCathode(
    redeclare package Medium = Medium_cath,
    T = T_cath_in,
    use_flow_in = true,
    use_Th_in = false,use_X_in = true) annotation (Placement(transformation(
          extent={{-102.0,-30.0},{-82.0,-10.0}},
                                      rotation=0.0,origin = {0.0,0.0})));
        
  .FuelCell.Sources.GasPressureBoundary sinkCathode(
    redeclare package Medium = Medium_cath,
    T = T_cath_in,use_p_in = true,use_Th_in = false,use_X_in = true) annotation (Placement(transformation(extent={{99.61496612607814,-30.38503387392182},{79.61496612607814,-10.38503387392182}},
          rotation=0.0,origin = {0.0,0.0})));
        
  .PEMFCModel4par.CoolStack coolStack(
    redeclare package Medium_an = Medium_an,
    redeclare package Medium_cath = Medium_cath,
    redeclare package Medium_cooling = MediumWater,
    X_feed_an = humidification_an.x_wet_in,
    X_feed_cath = humidification_cath.x_wet_in,p_start_out_anode = sinkAnode.p,p_start_in_anode = p_an_in,m_flow_start_cathode = flowCathode.m_flow,p_start_out_cathode = sinkCathode.p,p_start_in_cathode = p_cath_in,T_start_in_anode = flowAnode.T,T_start_in_cathode = flowCathode.T,T_start_out_cathode = flowCathode.T,T_start_out_anode = flowAnode.T,positiveFlow_cooling = true,from_dp_cooling = false,positiveFlow_anode = true,positiveFlow_cathode = true,from_dp_cathode = false,p_start_in_cooling = p_water_in,T_start_in_cooling = sourceW.T,T_start_out_cooling = sourceW.T,m_flow_start_cooling = sourceW.m_flow, h_inflow_an = flowAnode.fluidPort.h_outflow, h_inflow_cath = flowCathode.fluidPort.h_outflow, h_inflow_cooling = sourceW.fluidPort.h_outflow,p_start_out_cooling = sinkP.p,enable_setting = true,dp_smooth = 1e-2,mflow_smooth = 1e-8,dp_smooth_cooling = 1e-2,mflow_smooth_cooling = 1e-8,redeclare replaceable model Friction_anode = .Modelon.ThermoFluid.FlowChannels.PipeResistances.SinglePhase.LinearOperatingPointLoss(d0 = rho_an_nom,dp0 = dp_an_nom,m_flow0 = m_flow_an_nom),redeclare replaceable model Friction_cathode = .Modelon.ThermoFluid.FlowChannels.PipeResistances.SinglePhase.LinearOperatingPointLoss(d0 = rho_cath_nom,dp0 = dp_cath_nom,m_flow0 = m_flow_cath_nom),redeclare replaceable model HeatTransfer_anode = .Modelon.ThermoFluid.FlowChannels.HeatTransfer.SinglePhase.SinglePhase,redeclare replaceable model HeatTransfer_cathode = .Modelon.ThermoFluid.FlowChannels.HeatTransfer.SinglePhase.SinglePhase,redeclare replaceable model Friction_cooling = .Modelon.ThermoFluid.FlowChannels.PipeResistances.TwoPhase.SimpleFromDensity(dp0 = dp_water_nom,mflow0 = m_flow_water_nom),redeclare replaceable model HeatTransfer_cooling = .Modelon.ThermoFluid.FlowChannels.HeatTransfer.TwoPhase.SimpleTwoPhase,N = 5,A_cell = 300e-4,L_anode = 0.0813,D_anode = 0.023,D_cathode = 0.0445,L_cathode = 0.0813,L_cooling = 8.13,D_cooling = 0.0375,diffusiveSpecies = {"O2","N2"},n_cell = 11,M_stack = 10.1,z = 50e-6,j_0 = PID_j0.y,m_conc = PID_m_conc.y,alpha = 0.5,c1 = PID_c1.y,c_stack = 2e3,m_flow_start_anode = flowAnode.m_flow,X_start_anode = flowAnode.X,X_start_cathode = flowCathode.X,initOpt_anode = .Modelon.ThermoFluid.Choices.InitOptions.initialValues,initOpt_cooling = .Modelon.ThermoFluid.Choices.InitOptions.initialValues,initOpt_cathode = .Modelon.ThermoFluid.Choices.InitOptions.initialValues,n_channels_cooling = 1,from_dp_anode = false,CF_CathodeSidePressureLoss = 0.4,j_loss = PID_j_loss.y,n_conc = PID_n_conc.y) annotation (Placement(transformation(extent={{-24.0,-2.0},{16.0,38.0}},rotation = 0.0,origin = {0.0,0.0})));

   .FuelCell.Sources.WaterFlowBoundary sourceW(redeclare package Medium = MediumWater, T = T_water_in,use_flow_in = true,use_Th_in = false) annotation (Placement(transformation(extent={{-58.0,-50.0},{-38.0,-30.0}},rotation = 0.0,origin = {0.0,0.0})));
        
  .FuelCell.Sources.WaterPressureBoundary sinkP(redeclare package Medium = MediumWater, T = T_water_in,use_p_in = true) annotation (Placement(transformation(extent={{-10.0,10.0},{10.0,-10.0}},
        rotation=180.0,
        origin={48.58183823418737,-40.94544117720845})));
        
  .Modelon.Visualizers.RealValue display_V(number = coolStack.pin_p.v - coolStack.pin_n.v, precision = 3) annotation (Placement(transformation(extent={{-90,-86},{-70,-66}})));
  .Modelon.Visualizers.RealValue display_P(number = coolStack.pin_n.i*(coolStack.pin_p.v - coolStack.pin_n.v), precision = 0) annotation (Placement(transformation(extent={{-66,-86},{-46,-66}})));
  .Modelon.Visualizers.RealValue display_I(number = coolStack.pin_n.i, precision = 2) annotation (Placement(transformation(extent={{-90,-102},{-70,-82}})));
  .Modelon.Visualizers.RealValue display_j(number = coolStack.summary.j_external*1e-4*1e3, precision = 1) annotation (Placement(transformation(extent={{-38,-86},{-18,-66}})));
  .Modelon.Visualizers.RealValue display_Q(number = sum(coolStack.coolingPipe.q.Q_flow), precision = 1) annotation (Placement(transformation(extent={{-38,-102},{-18,-82}})));
  .Modelon.Visualizers.RealValue display_T(precision = 1, number = .Modelica.Units.Conversions.to_degC(coolStack.subStack.cell.T_cell_avg)) annotation (Placement(transformation(extent={{-68.0,-102.0},{-48.0,-82.0}},rotation = 0.0,origin = {0.0,0.0})));
  
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
    offset = m_flow_water,
    duration = 0,
    height = 0) annotation(Placement(transformation(extent = {{-82.70681533891384,-48.461718351330504},{-70.87211677910473,-36.62701979152139}},origin = {0.0,0.0},rotation = 0.0)));
  .Modelica.Electrical.Analog.Sources.ConstantCurrent constant_load(I = I) annotation(Placement(transformation(origin = {2.0,54.0},extent = {{7.0,-7.0},{-7.0,7.0}},rotation = -180.0)));
  .Modelica.Blocks.Sources.RealExpression dp_water(y = dp_water_nom) annotation(Placement(transformation(extent = {{20,-124},{40,-104}},origin = {0,0},rotation = 0)));
  .Modelica.Blocks.Math.Add p_an_out(k2 = +1,k1 = -1) annotation(Placement(transformation(extent = {{-52.0,118.0},{-32.0,138.0}},origin = {0.0,0.0},rotation = 0.0)));
  .Modelica.Blocks.Math.Add p_cath_out(k1 = +1,k2 = -1) annotation(Placement(transformation(extent = {{-54.0,-140.0},{-34.0,-120.0}},origin = {0.0,0.0},rotation = 0.0)));
  .Modelica.Blocks.Math.Add p_water_out(k2 = -1,k1 = +1) annotation(Placement(transformation(extent = {{62.0,-104.0},{82.0,-84.0}},origin = {0.0,0.0},rotation = 0.0)));
    .Modelica.Electrical.Analog.Sensors.VoltageSensor voltageSensor annotation(Placement(transformation(extent = {{-1.966451765551266,119.96645176555127},{9.966451765551266,108.03354823444873}},origin = {0.0,0.0},rotation = 0.0)));
    .Modelica.Blocks.Math.Gain gain(k = 1 / coolStack.n_cell) annotation(Placement(transformation(extent = {{-606.344731741886,53.655268258114035},{-593.655268258114,66.34473174188597}},origin = {0.0,0.0},rotation = 0.0)));
    .Modelica.Blocks.Continuous.LimPID PID_j0(controllerType = .Modelica.Blocks.Types.SimpleController.P,yMin = 1,yMax = 10,k = 1e-1,initType = .Modelica.Blocks.Types.Init.InitialState) annotation(Placement(transformation(extent = {{-524.0,26.0},{-504.0,46.0}},origin = {0.0,0.0},rotation = 0.0)));
    .Modelica.Blocks.Continuous.LimPID PID_j_loss(controllerType = .Modelica.Blocks.Types.SimpleController.P,yMin = 1,yMax = 10,k = 1e-1,initType = .Modelica.Blocks.Types.Init.InitialState) annotation(Placement(transformation(extent = {{-524.0,-32.0},{-504.0,-12.0}},origin = {0.0,0.0},rotation = 0.0)));
    .Modelica.Blocks.Continuous.LimPID PID_m_conc(controllerType = .Modelica.Blocks.Types.SimpleController.P,yMin = 1e-6,yMax = 1,k = 1e-1,initType = .Modelica.Blocks.Types.Init.InitialState) annotation(Placement(transformation(extent = {{-522.0,-100.0},{-502.0,-80.0}},origin = {0.0,0.0},rotation = 0.0)));
    .Modelica.Blocks.Continuous.LimPID PID_n_conc(controllerType = .Modelica.Blocks.Types.SimpleController.P,yMin = 1e-6,yMax = 1,k = 1e-1,initType = .Modelica.Blocks.Types.Init.InitialState) annotation(Placement(transformation(extent = {{-522.0,-166.0},{-502.0,-146.0}},origin = {0.0,0.0},rotation = 0.0)));
    .Modelica.Blocks.Interaction.Show.RealValue relative_error(significantDigits = 3,use_numberPort = false,number = err_rel) annotation(Placement(transformation(extent = {{-24.08338897633281,-208.0833889763328},{12.08338897633281,-171.9166110236672}},origin = {0.0,0.0},rotation = 0.0)));
    .Modelica.Blocks.Interaction.Show.RealValue V_cell(number = coolStack.subStack.cell.V_cell,use_numberPort = false,significantDigits = 3) annotation(Placement(transformation(extent = {{-24.083388976332813,-238.0833889763328},{12.083388976332813,-201.9166110236672}},origin = {0.0,0.0},rotation = 0.0)));
    .Modelica.Blocks.Sources.RealExpression dp_an_SS(y = dp_an) annotation(Placement(transformation(extent = {{-144.0,138.0},{-124.0,158.0}},origin = {0.0,0.0},rotation = 0.0)));
    .Modelica.Blocks.Sources.RealExpression p_an_SS(y = p_an_in) annotation(Placement(transformation(extent = {{-200.0,68.0},{-180.0,88.0}},origin = {0.0,0.0},rotation = 0.0)));
    .Modelica.Blocks.Sources.RealExpression flow_rate_an_SS(y = V_flow_an) annotation(Placement(transformation(extent = {{-200.0,28.0},{-180.0,48.0}},origin = {0.0,0.0},rotation = 0.0)));
    .Modelica.Blocks.Sources.RealExpression flow_rate_cath_SS(y = V_flow_cath) annotation(Placement(transformation(extent = {{-196.0,-44.0},{-176.0,-24.0}},origin = {0.0,0.0},rotation = 0.0)));
    .Modelica.Blocks.Sources.RealExpression p_cath_SS(y = p_cath_in) annotation(Placement(transformation(extent = {{-200.0,-82.0},{-180.0,-62.0}},origin = {0.0,0.0},rotation = 0.0)));
    .Modelica.Blocks.Sources.RealExpression dp_cath_SS(y = dp_cath) annotation(Placement(transformation(extent = {{-188.0,-166.0},{-168.0,-146.0}},origin = {0.0,0.0},rotation = 0.0)));
    .Modelica.Blocks.Sources.RealExpression p_water_SS(y = p_water_in) annotation(Placement(transformation(extent = {{20.0,-90.0},{40.0,-70.0}},origin = {0.0,0.0},rotation = 0.0)));
    .Modelica.Blocks.Sources.RealExpression realExpression(y = V_exp) annotation(Placement(transformation(extent = {{-608.0,92.0},{-588.0,112.0}},origin = {0.0,0.0},rotation = 0.0)));
    .Modelica.Blocks.Continuous.LimPID PID_c1(controllerType = .Modelica.Blocks.Types.SimpleController.P,yMax = 1e-2,yMin = 1e-6,initType = .Modelica.Blocks.Types.Init.InitialState,k = 1e-2) annotation(Placement(transformation(extent = {{-524,92},{-504,112}},origin = {0,0},rotation = 0)));
           
equation       

  connect(multiDisplaySensor4.u, display_phTmdot_w_in. y) annotation (Line(points={{-25.95,-39.95},{-25.547635814819174,-32}}, color={0,0,0}));
  connect(multiDisplaySensor1.u, display_phTmdot_w_out. y) annotation (Line(points={{24.05,
          -39.95},{24.05,-35.975},{24,-35.975},{24,-32}}, color={0,0,0}));
  connect(multiDisplaySensor4.portA, sourceW.fluidPort) annotation (Line(points={{-32,-40},{-39,-40}}, color={0,0,255}));
  connect(multiDisplaySensor4.portB, coolStack.feed_cooling) annotation (Line(points={{-20,-40},{-12,-40},{-12,1.5}}, color={0,0,255}));
  connect(multiDisplaySensor1.portB, sinkP.fluidPort) annotation (Line(points={{30,-40},{39.58183823418737,-40},{39.58183823418737,-40.94544117720845}}, color={0,0,255}));
  connect(multiDisplaySensor1.portA, coolStack.drain_cooling) annotation (Line(points={{18,-40},{4,-40},{4,1.5}}, color={0,0,255}));
  connect(display_phTmdot_an_in.y, gasSensor3. u) annotation (Line(points={{-68,40},{
          -68,30.05},{-67.95,30.05}},
                                  color={0,0,0}));
  connect(display_an_in.data, gasSensor3. u) annotation (Line(points={{-40.462753249644805,37.54955698575693},{-56,37.54955698575693},{-56,30.05},{-67.95,30.05}},
                                                         color={0,0,0}));
  connect(gasSensor3.portA, flowAnode.fluidPort) annotation (Line(points={{-74,30},{-82.93906632502666,30},{-82.93906632502666,29.768425703538625}}, color={255,128,0}));
  connect(gasSensor3.portB, coolStack.feed_an) annotation (Line(points={{-62,30},{-22,30},{-22,28}}, color={255,128,0}));
  connect(gasSensor2.u, display_phTmdot_an_out. y) annotation (Line(points={{34.05,
          30.05},{34.05,35.025},{34,35.025},{34,40}},
                                               color={0,0,0}));
  connect(gasSensor2.u, display_an_out.data) annotation (Line(points={{34.05,30.05},{34.05,32},{62,32},{62,37.2941}},
                                                   color={0,0,0}));
  connect(gasSensor2.portA, coolStack.drain_an) annotation (Line(points={{28,30},{14,30},{14,28}}, color={255,128,0}));
  connect(gasSensor2.portB, sinkAnode.fluidPort) annotation (Line(points={{40,30},{80.22993225215629,30},{80.22993225215629,29.6149661260782}}, color={255,128,0}));
  connect(gasSensor4.u, display_phTmdot_cath_in. y) annotation (Line(points={{-67.95,-19.95},{-67.95,-14.975},{-68,-14.975},{-68,-10}},
                                                     color={0,0,0}));
  connect(gasSensor4.u, display_cath_in. data) annotation (Line(points={{-67.95,-19.95},{-42.744560661301904,-19.95},{-42.744560661301904,-14}}, color={0,0,0}));
  connect(gasSensor4.portA, flowCathode.fluidPort) annotation (Line(points={{-74,-20},{-83,-20}}, color={255,128,0}));
  connect(gasSensor4.portB, coolStack.feed_cath) annotation (Line(points={{-62,-20},{-22,-20},{-22,8}}, color={255,128,0}));
  connect(gasSensor1.u, display_cath_out.data) annotation (Line(points={{36.05,-19.95},{62.00000000000004,-19.95},{62.00000000000004,-13}},    color={0,0,0}));
  connect(sinkCathode.fluidPort, gasSensor1.portB) annotation (Line(points={{80.61496612607814,-20.38503387392182},{80.61496612607814,-20},{42,-20}}, color={255,128,0}));
  connect(gasSensor1.portA, coolStack.drain_cath) annotation (Line(points={{30,-20},{14,-20},{14,8}}, color={255,128,0}));
  connect(gasSensor1.u, display_phTmdot_cath_out.y) annotation (Line(points={{36.05,
          -19.95},{36.05,-14.975},{36,-14.975},{36,-10}}, color={0,0,0}));
  connect(humidification_an.m_flow_wet_in,flowAnode.m_flow_in) annotation(Line(points = {{-116.40952305789519,32.588021244293536},{-110.40952305789519,32.588021244293536},{-110.40952305789519,45.768425703538625},{-97.93906632502666,45.768425703538625},{-97.93906632502666,39.768425703538625}},color = {0,0,127}));
  connect(humidification_an.x_wet_in,flowAnode.X_in) annotation(Line(points = {{-116.18859441508964,27.91695851069057},{-110.18859441508964,27.91695851069057},{-110.18859441508964,45.768425703538625},{-85.93906632502666,45.768425703538625},{-85.93906632502666,39.768425703538625}},color = {0,0,127}));
  connect(humidification_cath.m_flow_wet_in,flowCathode.m_flow_in) annotation(Line(points = {{-116.40952305789519,-15.411978755706466},{-110.40952305789519,-15.411978755706466},{-110.40952305789519,-4},{-98,-4},{-98,-10}},color = {0,0,127}));
  connect(humidification_cath.x_wet_in,flowCathode.X_in) annotation(Line(points = {{-116.18859441508964,-20.08304148930943},{-110,-20.08304148930943},{-110,-4},{-86,-4},{-86,-10}},color = {0,0,127}));
  connect(ramp_m_flow_water.y,sourceW.m_flow_in) annotation(Line(points = {{-70.28038185111427,-42.54436907142595},{-61.95890105565714,-42.54436907142595},{-61.95890105565714,-30},{-54,-30}},color = {0,0,127}));
  connect(humidification_an.x_wet_in,sinkAnode.X_in) annotation(Line(points = {{-116.18859441508964,27.91695851069057},{-110.18859441508964,27.91695851069057},{-110.18859441508964,84},{83.22993225215629,84},{83.22993225215629,39.6149661260782}},color = {0,0,127}));
  connect(humidification_cath.x_wet_in,sinkCathode.X_in) annotation(Line(points = {{-116.18859441508964,-20.08304148930943},{-110,-20.08304148930943},{-110,16},{83.61496612607814,16},{83.61496612607814,-10.38503387392182}},color = {0,0,127}));
  connect(coolStack.pin_n,ground.p) annotation(Line(points = {{4,34.5},{16,34.5},{16,64},{30,64}},color = {0,0,255}));      
  connect(constant_load.p,coolStack.pin_p) annotation(Line(points = {{-5.000000000000001,54},{-10,54},{-10,34.5}},color = {0,0,255}));
  connect(constant_load.n,coolStack.pin_n) annotation(Line(points = {{9,54},{16,54},{16,34.5},{6,34.5}},color = {0,0,255})); 
  connect(p_an_out.y,sinkAnode.p_in) annotation(Line(points = {{-31,128},{95.22993225215629,128},{95.22993225215629,39.6149661260782}},color = {0,0,127}));
  connect(p_cath_out.y,sinkCathode.p_in) annotation(Line(points = {{-33,-130},{112,-130},{112,-4.38503387392182},{95.61496612607814,-4.38503387392182},{95.61496612607814,-10.38503387392182}},color = {0,0,127}));
  connect(p_water_out.y,sinkP.p_in) annotation(Line(points = {{83,-94},{92,-94},{92,-26},{54.58183823418737,-26},{54.58183823418737,-30.94544117720845}},color = {0,0,127}));
  connect(dp_water.y,p_water_out.u2) annotation(Line(points = {{41,-114},{50.5,-114},{50.5,-100},{60,-100}},color = {0,0,127}));
    connect(voltageSensor.p,coolStack.pin_p) annotation(Line(points = {{-1.966451765551266,114},{-12,114},{-12,34.5}},color = {0,0,255}));
    connect(voltageSensor.n,coolStack.pin_n) annotation(Line(points = {{9.966451765551266,114},{16,114},{16,34.5},{4,34.5}},color = {0,0,255}));
    connect(voltageSensor.v,gain.u) annotation(Line(points = {{4,120.5630969421064},{4,198},{-658,198},{-658,60},{-607.6136780902632,60}},color = {0,0,127}));
    connect(gain.y,PID_j0.u_m) annotation(Line(points = {{-593.0207950839253,60},{-589.655268258114,60},{-589.655268258114,18},{-514,18},{-514,24}},color = {0,0,127}));
    connect(gain.y,PID_j_loss.u_m) annotation(Line(points = {{-593.0207950839253,60},{-590,60},{-590,-40},{-514,-40},{-514,-34}},color = {0,0,127}));
    connect(gain.y,PID_m_conc.u_m) annotation(Line(points = {{-593.0207950839253,60},{-590,60},{-590,-108},{-512,-108},{-512,-102}},color = {0,0,127}));
    connect(PID_n_conc.u_m,gain.y) annotation(Line(points = {{-512,-168},{-512,-174},{-590,-174},{-590,60},{-593.02,60}},color = {0,0,127}));
    connect(dp_an_SS.y,p_an_out.u1) annotation(Line(points = {{-123,148},{-88.5,148},{-88.5,134},{-54,134}},color = {0,0,127}));
    connect(p_an_SS.y,p_an_out.u2) annotation(Line(points = {{-179,78},{-146,78},{-146,122},{-54,122}},color = {0,0,127}));
    connect(p_an_SS.y,humidification_an.p_in) annotation(Line(points = {{-179,78},{-146,78},{-146,33.78734816238078},{-131.57469632476156,33.78734816238078}},color = {0,0,127}));
    connect(flow_rate_an_SS.y,humidification_an.V_flow_dry_in) annotation(Line(points = {{-179,38},{-155.2873481623808,38},{-155.2873481623808,26.21265183761922},{-131.57469632476156,26.21265183761922}},color = {0,0,127}));
    connect(flow_rate_cath_SS.y,humidification_cath.V_flow_dry_in) annotation(Line(points = {{-175,-34},{-157.2873481623808,-34},{-157.2873481623808,-21.78734816238078},{-131.57469632476156,-21.78734816238078}},color = {0,0,127}));
    connect(p_cath_SS.y,humidification_cath.p_in) annotation(Line(points = {{-179,-72},{-146,-72},{-146,-14.212651837619218},{-131.57469632476156,-14.212651837619218}},color = {0,0,127}));
    connect(p_cath_SS.y,p_cath_out.u1) annotation(Line(points = {{-179,-72},{-146,-72},{-146,-124},{-56,-124}},color = {0,0,127}));
    connect(dp_cath_SS.y,p_cath_out.u2) annotation(Line(points = {{-167,-156},{-111.5,-156},{-111.5,-136},{-56,-136}},color = {0,0,127}));
    connect(p_water_SS.y,p_water_out.u1) annotation(Line(points = {{41,-80},{50.5,-80},{50.5,-88},{60,-88}},color = {0,0,127}));
    connect(realExpression.y,PID_j0.u_s) annotation(Line(points = {{-587,102},{-556.5,102},{-556.5,36},{-526,36}},color = {0,0,127}));
    connect(realExpression.y,PID_j_loss.u_s) annotation(Line(points = {{-587,102},{-556,102},{-556,-22},{-526,-22}},color = {0,0,127}));
    connect(realExpression.y,PID_m_conc.u_s) annotation(Line(points = {{-587,102},{-556,102},{-556,-90},{-524,-90}},color = {0,0,127}));
    connect(realExpression.y,PID_n_conc.u_s) annotation(Line(points = {{-587,102},{-556,102},{-556,-156},{-524,-156}},color = {0,0,127}));
    connect(gain.y,PID_c1.u_m) annotation(Line(points = {{-593.0207950839253,60},{-514,60},{-514,90}},color = {0,0,127}));
    connect(realExpression.y,PID_c1.u_s) annotation(Line(points = {{-587,102},{-526,102}},color = {0,0,127}));
    
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

end CoolStack_PID5_constantcurrent;
