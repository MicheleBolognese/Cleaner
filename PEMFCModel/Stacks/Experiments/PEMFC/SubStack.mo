within PEMFCModel.Stacks.Experiments.PEMFC;
model SubStack "PEMFC substack test model"
    
  extends .Modelon.Icons.Experiment;
    
  replaceable package Medium_an = .PEMFCModel4par.Medium.NASASteamH2N2O2 constrainedby
    .FuelCell.Media.Templates.ReactionGas "Anodic medium";
  replaceable package Medium_cath = .FuelCell.Media.PreDefined.IdealGases.NASAMoistAir constrainedby
    .FuelCell.Media.Templates.ReactionGas "Cathodic medium";

  parameter .Modelica.Units.SI.MassFlowRate m_flow_an  = 2.901e-4 "Anodic gas mass flow rate" annotation(Dialog(group = "Anode side"));
  parameter .Modelica.Units.SI.Temperature T_an = 70 + 273.15 "Anodic gas temperature" annotation(Dialog(group = "Anode side"));
  parameter .Modelica.Units.SI.Pressure p_an_in_start = 1.9e5 "Start value of inlet gas pressure at anode" annotation(Dialog(group = "Anode side"));
  parameter .Modelica.Units.SI.Pressure p_an_out = 1.8e5 "Anodic gas outlet pressure" annotation(Dialog(group = "Anode side"));
  parameter .Modelica.Units.SI.MoleFraction y_an[Medium_an.nS] = {0.589,0.158,0.253,0} "Anodic gas molar-based composition" annotation(Dialog(group = "Anode side"));
  
  parameter .Modelica.Units.SI.MassFlowRate m_flow_cath = 1.096e-3 "Cathodic gas mass flow rate" annotation(Dialog(group = "Cathode side"));
  parameter .Modelica.Units.SI.Temperature T_cath = 70 + 273.15 "Cathodic gas temeprature" annotation(Dialog(group = "Cathode side"));
  parameter .Modelica.Units.SI.Pressure p_cath_in_start = 1.7e5 "Start value of inlet gas pressure at cathode" annotation(Dialog(group = "Cathode side"));
  parameter .Modelica.Units.SI.Pressure p_cath_out = 1.55e5 "Cathodic gas outlet pressure" annotation(Dialog(group = "Cathode side"));
  parameter .Modelica.Units.SI.MoleFraction y_cath[Medium_cath.nS] = {0,0,0.204,0.629,0.167}  "Cathodic gas molar-based composition" annotation(Dialog(group = "Cathode side"));
  
  parameter .Modelica.Units.SI.HeatFlowRate Q_cool_min = -139 "Minimum value of cooling power" annotation(Dialog(group = "Cooling"));
  parameter .Modelica.Units.SI.HeatFlowRate Q_cool_max = Q_cool_min "Maximum value of cooling power" annotation(Dialog(group = "Cooling"));
    
  .FuelCell.Sources.GasFlowBoundary flowCathode(
    redeclare package Medium = Medium_cath,
    m_flow = m_flow_cath,
    T = T_cath,
    X = Medium_cath.moleToMassFractions(y_cath,Medium_cath.MMX)) annotation (Placement(transformation(
          extent={{-100.0,-20.0},{-80.0,0.0}},rotation=0.0,origin = {0.0,0.0})));
  
  .FuelCell.Sources.GasFlowBoundary flowAnode(
    redeclare package Medium = Medium_an,
    m_flow = m_flow_an,
    X = Medium_an.moleToMassFractions(y_an,Medium_an.MMX),
    T = T_an) annotation (Placement(transformation(
          extent={{-100.0,17.774919282245833},{-80.0,37.77491928224583}},
                                     rotation=0.0,origin = {0.0,0.0})));

  .FuelCell.Sources.GasPressureBoundary sinkCathode(
    redeclare package Medium = Medium_cath,
    p = p_cath_out,
    T = T_cath,X = Medium_cath.moleToMassFractions(y_cath,Medium_cath.MMX)) annotation (Placement(transformation(extent={{98.0,-20.0},{78.0,0.0}},
          rotation=0.0,origin = {0.0,0.0})));
    
  .FuelCell.Sources.GasPressureBoundary sinkAnode(
    redeclare package Medium = Medium_an,
    p = p_an_out,
    T = T_an,X = Medium_an.moleToMassFractions(y_an,Medium_an.MMX)) annotation (Placement(transformation(extent={{98,18},{78,38}},
          rotation=0)));

  PEMFCModel.Stacks.PEMFC.Empirical.SubStack subStack(
    redeclare package Medium_an = Medium_an,
    redeclare package Medium_cath = Medium_cath,
    N = 5,
    n_cell = 11,
    h_inflow_an = flowAnode.fluidPort.h_outflow,
    h_inflow_cath = flowCathode.fluidPort.h_outflow,
    X_feed_an = Medium_an.moleToMassFractions(y_an,Medium_an.MMX),
    X_feed_cath = Medium_cath.moleToMassFractions(y_cath,Medium_cath.MMX),
    positiveFlow_anode = true,positiveFlow_cathode = true,T_start_in_anode = flowAnode.T,T_start_out_anode = flowAnode.T,X_start_anode = flowAnode.X,m_flow_start_anode = m_flow_an,p_start_out_cathode = sinkCathode.p,p_start_in_anode = p_an_in_start,p_start_out_anode = sinkAnode.p,p_start_in_cathode = p_cath_in_start,T_start_in_cathode = flowCathode.T,T_start_out_cathode = flowCathode.T,X_start_cathode = flowCathode.X,m_flow_start_cathode = m_flow_cath,from_dp_anode = false,from_dp_cathode = false,redeclare replaceable model HeatTransfer_anode = .Modelon.ThermoFluid.FlowChannels.HeatTransfer.SinglePhase.SinglePhase,redeclare replaceable model HeatTransfer_cathode = .ThermalPower.SubComponents.HeatTransfer.SinglePhase.Plate.SinglePhase,redeclare replaceable model Friction_anode = .Modelon.ThermoFluid.FlowChannels.PipeResistances.SinglePhase.LinearOperatingPointLoss(dp0 = 0.1e5,m_flow0 = m_flow_an,d0 = 0.7),redeclare replaceable model Friction_cathode = .Modelon.ThermoFluid.FlowChannels.PipeResistances.SinglePhase.LinearOperatingPointLoss(d0 = 1.6,dp0 = 0.15e5,m_flow0 = m_flow_cath),initFromEnthalpy_anode = false,mflow_smooth = 1e-8,dp_smooth = 1e-2,A_cell = 300e-4,M_stack = 10.1,L_anode = 0.0813,D_anode = 0.023,L_cathode = 0.0813,D_cathode = 0.0445,z = 75e-6,diffusiveSpecies = {"O2","N2"},m_conc = 0,j_loss = 0) annotation (Placement(transformation(extent={{-17.458305704687554,1.6250828859373563},{22.541694295312446,41.62508288593736}}, rotation=0.0,origin = {0.0,0.0})));

  .Modelica.Electrical.Analog.Basic.Ground ground annotation (Placement(transformation(extent={{-32.60355701537692,-39.99435108118894},{-16.603557015376918,-23.994351081188942}},   rotation=
           0.0,origin = {0.0,0.0})));
    
  .Modelica.Electrical.Analog.Sources.RampCurrent current(
    startTime = 3600,
    offset = 180,
    I = 0,
    duration = 0) annotation (Placement(transformation(
        origin={0.0,0.0},
        extent={{2.0,-28.0},{-18.0,-8.0}},
        rotation=0.0)));
    
  .Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow cooling[subStack.N](each alpha = 0, each T_ref = .Modelica.Units.Conversions.from_degC(800)) annotation (Placement(transformation(
        origin={8.0,52.0},
        extent={{-10.0,-10.0},{10.0,10.0}},
        rotation=-90.0)));
  .Modelica.Blocks.Sources.Ramp Q_flow[subStack.N](
    each startTime = 3600,
    each duration = 0,
    each height = Q_cool_max - Q_cool_min,offset = fill(Q_cool_min,subStack.N)) annotation (Placement(transformation(extent={{-18.0,64.0},{-2.0,80.0}},         rotation=0.0,origin = {0.0,0.0})));
 
  .Modelon.Visualizers.RealValue display_V(number = current.p.v - current.n.v, precision = 3) annotation (Placement(transformation(extent={{-90,-86},{-70,-66}})));
  .Modelon.Visualizers.RealValue display_P(number = current.p.i*(current.p.v - current.n.v), precision = 0) annotation (Placement(transformation(extent={{-66,-86},{-46,-66}})));
  .Modelon.Visualizers.RealValue display_I(number = current.p.i, precision = 2) annotation (Placement(transformation(extent={{-90,-102},{-70,-82}})));
  .Modelon.Visualizers.RealValue display_j(precision = 1, number = subStack.summary.j_external*1e-4*1e3) annotation (Placement(transformation(extent={{-38,-86},{-18,-66}})));
  .Modelon.Visualizers.RealValue display_Q(precision = 1, number = -sum(subStack.wall.Q_flow)) annotation (Placement(transformation(extent={{-38,-102},{-18,-82}})));
  .Modelon.Visualizers.RealValue display_T(precision = 1, number = .Modelica.Units.Conversions.to_degC(subStack.cell.T_cell_avg)) annotation (Placement(transformation(extent={{-66,-102},{-46,-82}})));
  
  .FuelCell.Sensors.GasMultiDisplaySensor gasSensor3(redeclare package Medium = Medium_an) annotation (Placement(transformation(extent={{-76,18},{-56,38}})));
  .FuelCell.Sensors.GasMultiDisplaySensor gasSensor4(redeclare package Medium = Medium_cath) annotation (Placement(transformation(extent={{-76,-20},{-56,0}})));
  .FuelCell.Sensors.GasMultiDisplaySensor gasSensor2(redeclare package Medium = Medium_an) annotation (Placement(transformation(extent={{28.29027197116028,18.0},{48.29027197116028,38.0}},rotation = 0.0,origin = {0.0,0.0})));
  .FuelCell.Sensors.GasMultiDisplaySensor gasSensor1(redeclare package Medium = Medium_cath) annotation (Placement(transformation(extent={{32,-20},{52,0}})));
 
  .FuelCell.Sensors.MultiDisplay_phTmdot display_phTmdot_an_in( displayUnits = true, redeclare package Medium = Medium_an) annotation (Placement(transformation(extent={{-76.76631801609432,28.0},{-56.76631801609432,48.0}},rotation = 0.0,origin = {0.0,0.0})));
  .FuelCell.Sensors.MultiDisplay_phTmdot display_phTmdot_cath_in(displayUnits = true, redeclare package Medium = Medium_cath) annotation (Placement(transformation(extent={{-76,-10},{-56,10}})));
  .FuelCell.Sensors.MultiDisplay_phTmdot display_phTmdot_an_out(displayUnits = true, redeclare package Medium = Medium_an) annotation (Placement(transformation(extent={{30.0,30.0},{50.0,50.0}},rotation = 0.0,origin = {0.0,0.0})));
  .FuelCell.Sensors.MultiDisplay_phTmdot display_phTmdot_cath_out(displayUnits = true, redeclare package Medium = Medium_cath) annotation (Placement(transformation(extent={{32,-12},{52,8}})));
    
  .FuelCell.Sensors.ReformateLongCompositionDisplay display_MoleFractions_an_in(redeclare package Medium = Medium_an) annotation (Placement(transformation(extent={{-51.4891213226038,40.766318016094324},{-27.4891213226038,60.766318016094324}},rotation = 0.0,origin = {0.0,0.0})));
  .FuelCell.Sensors.AirCompositionDisplay display_MoleFractions_cath_in(redeclare package Medium = Medium_cath) annotation (Placement(transformation(extent={{-52,-4},{-30,14}})));
  .FuelCell.Sensors.ReformateLongCompositionDisplay display_MoleFractions_an_out(redeclare package Medium = Medium_an) annotation (Placement(transformation(extent={{52,38},{76,58}})));
  .FuelCell.Sensors.AirCompositionDisplay display_MoleFractions_cath_out(redeclare package Medium = Medium_cath) annotation (Placement(transformation(extent={{56,-6},{76,12}})));
   
 
  
    
equation
  
  connect(current.n, ground. p) annotation (Line(
      points={{-18,-18},{-24.603557015376918,-18},{-24.603557015376918,-23.994351081188942}},
      color={0,0,255},
      smooth=Smooth.None));
  connect(current.n, subStack.pin_n) annotation (Line(
      points={{-18,-18},{-24.80418380686819,-18},{-24.80418380686819,39.443762417327775},{2.541694295312446,39.443762417327775},{2.541694295312446,34.825082885937356}},
      color={0,0,255},
      smooth=Smooth.None));
  connect(current.p, subStack.pin_p) annotation (Line(points={{2,-18},{2,8.425082885937355},{2.541694295312446,8.425082885937355}},               color={0,0,255}));
  
  connect(Q_flow.y, cooling.Q_flow) annotation (Line(
      points={{-1.1999999999999993,72},{7.999999999999998,72},{7.999999999999998,62}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(cooling.port, subStack.wall) annotation (Line(
      points={{8.000000000000002,42},{8.000000000000002,21.625082885937356},{10.541694295312446,21.625082885937356}},
      color={191,0,0},
      smooth=Smooth.None));
  
  connect(display_phTmdot_an_in.y,gasSensor3. u) annotation (Line(points={{-66.76631801609432,38},{-66.76631801609432,28.05},{-65.95,28.05}},
                                  color={0,0,0}));
  connect(display_MoleFractions_an_in.data,gasSensor3. u) annotation (Line(points={{-39.4891213226038,36.06043566315315},{-39.4891213226038,28.05},{-65.95,28.05}}, color={0,0,0}));
  connect(gasSensor3.portA, flowAnode.fluidPort) annotation (Line(points={{-72,28},{-81,28},{-81,27.774919282245833}}, color={255,128,0}));
  connect(gasSensor3.portB, subStack.feed_anode) annotation (Line(points={{-60,28},{-15.458305704687554,28},{-15.458305704687554,29.625082885937356}}, color={255,128,0}));
  connect(gasSensor4.u,display_phTmdot_cath_in. y) annotation (Line(points={{-65.95,
          -9.95},{-65.95,-4.975},{-66,-4.975},{-66,0}},
                                                     color={0,0,0}));
  connect(gasSensor4.u,display_MoleFractions_cath_in. data) annotation (Line(points={{-65.95,
          -9.95},{-65.95,-7.975},{-41,-7.975},{-41,-4}},     color={0,0,0}));
  connect(subStack.feed_cathode, gasSensor4.portB) annotation (Line(points={{-15.458305704687554,13.625082885937356},{-15.458305704687554,-10},{-60,-10}}, color={255,128,0}));
  connect(gasSensor4.portA, flowCathode.fluidPort) annotation (Line(points={{-72,-10},{-81,-10}}, color={255,128,0}));
  connect(display_MoleFractions_an_out.data, gasSensor2.u) annotation (Line(points={{64,33.2941},{64,28.05},{38.34027197116028,28.05}},    color={0,0,0}));
  connect(gasSensor2.u, display_phTmdot_an_out.y) annotation (Line(points={{38.34027197116028,28.05},{38.34027197116028,32.025},{40,32.025},{40,40}}, color={0,0,0}));
  connect(sinkAnode.fluidPort, gasSensor2.portB) annotation (Line(points={{79,28},{44.29027197116028,28}}, color={255,128,0}));
  connect(gasSensor2.portA, subStack.drain_anode) annotation (Line(points={{32.29027197116028,28},{20.541694295312446,28},{20.541694295312446,29.625082885937356}}, color={255,128,0}));
  connect(gasSensor1.u,display_phTmdot_cath_out. y) annotation (Line(points={{42.05,
          -9.95},{42.05,-4.975},{42,-4.975},{42,-2}},
                                                  color={0,0,0}));
  connect(gasSensor1.u, display_MoleFractions_cath_out.data) annotation (Line(points={{
          42.05,-9.95},{42.05,-8},{66,-8},{66,-6}}, color={0,0,0}));
  connect(sinkCathode.fluidPort, gasSensor1.portB) annotation (Line(points={{79,-10},{48,-10}}, color={255,128,0}));
  connect(gasSensor1.portA, subStack.drain_cathode) annotation (Line(points={{36,-10},{20.541694295312446,-10},{20.541694295312446,13.625082885937356}}, color={255,128,0}));
    
  annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-120,
            -100},{100,100}}), graphics={
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
    experiment(StopTime=6000),
    __Dymola_experimentSetupOutput(equdistant=true),
    Documentation(revisions="<html>
Copyright &copy; 2004-2024, MODELON AB <br /> The use of this software component is regulated by the licensing conditions for Modelon Libraries. <br />This copyright notice must, unaltered, accompany all components that are derived from, copied from, <br />or by other means have their origin from any Modelon Library.
</html>", info="<html>
<h4>TestPEMFCSubstack</h4>
<p>Experiment with proton exchange membrane fuel cell substack model. The current is ramped from 10 to 60 ampere during 4000 seconds, which results in an increased heat production from the cell.
Note that the cooling flow is set to zero.</p>
</html>"),
    Icon(coordinateSystem(extent={{-120,-100},{100,100}})));
    
end SubStack;
