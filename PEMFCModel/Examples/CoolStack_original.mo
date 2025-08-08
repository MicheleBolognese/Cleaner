within PEMFCModel.Examples;

model CoolStack_original "Cooled PEMFC stack test model"
  extends .Modelon.Icons.Experiment;

  package Medium_an =
      .PEMFCModel4par.Medium.NASASteamH2N2O2
    "Medium at anode side";
  package Medium_cath = .FuelCell.Media.PreDefined.IdealGases.NASAMoistAir
    "Medium at cathode side";
  package MediumWater = .Modelon.Media.PreDefined.TwoPhase.WaterIF97
    "Water medium";

  parameter Integer[1] cath = Medium_cath.substanceIndexVector({"O2"})
    "index of required substances on the cathode side"                                                                        annotation(Evaluate=true);

  parameter .Modelica.Units.SI.MassFraction[Medium_an.nS] X_fuel=Medium_an.moleToMassFractions({0.589,0.158,0.253,0},Medium_an.MMX)
    "Inlet fuel composition";

  parameter .Modelica.Units.SI.MassFraction[Medium_cath.nS] X_air=Medium_cath.moleToMassFractions({0,0,0.204,0.629,0.167},Medium_cath.MMX) "Inlet air composition";

  parameter .Modelica.Units.SI.Temperature Tin_fuel=70 + 273.15 "Inlet fuel temperature";
  parameter .Modelica.Units.SI.Temperature Tin_air=70 + 273.15 "Inlet air temperature";
  parameter .Modelica.Units.SI.Temperature Tin_water=65 + 273.15 "Inlet water temperature";

  parameter .Modelica.Units.SI.Pressure pstart_an=1.9e5 "Start pressure (anode side)";
  parameter .Modelica.Units.SI.Pressure pstart_cath=1.7e5 "Start pressure (cathode side)";
  parameter .Modelica.Units.SI.Pressure pin_water=2.1e5 "Inlet water pressure";

  parameter .Modelica.Units.SI.MassFlowRate m_fuel=2.901e-4 "Inlet fuel flow rate";
  parameter .Modelica.Units.SI.MassFlowRate m_air=1.096e-3 "Inlet air flow rate";
  parameter .Modelica.Units.SI.MassFlowRate m_water=8.25e-2 "Inlet water flow rate";
  parameter Real cath_stoich = 1.801 "cathode stoichiometry";
  parameter Real anode_stoich = 1.5 "anode stoichiometry";

  parameter .Modelica.Units.SI.MassFlowRate m_out=m_air*X_air[cath[1]]*(1 - 1/cath_stoich)
    "calculated mass flow rate of oxygen gas at cathode outlet, may be used for control of the inlet cathode flow rate";
  .FuelCell.Sources.GasFlowBoundary                      flowAnode(
    redeclare replaceable package Medium = Medium_an,
    m_flow=m_fuel,
    X=X_fuel,
    T=Tin_fuel,
    use_flow_in=true,
    use_Th_in=false)                      annotation (Placement(transformation(
          extent={{-100,20},{-80,40}},
                                     rotation=0)));

  .FuelCell.Sources.GasPressureBoundary             sinkAnode(
    X=X_fuel,
    p=1.8e5,
    T=Tin_air,redeclare replaceable package Medium = Medium_an)
                   annotation (Placement(transformation(extent={{100,20},{80,40}},
          rotation=0)));

  .Modelica.Electrical.Analog.Basic.Ground ground
    annotation (Placement(transformation(extent={{-26,50},{-14,62}},   rotation=
           0)));
  .Modelica.Electrical.Analog.Sources.RampCurrent current(
    startTime=0,
    I=0,
    duration=3600,
    offset=180)
    annotation (Placement(transformation(
        origin={1.0,69.0},
        extent={{-7.0,-7.0},{7.0,7.0}},
        rotation=180.0)));
  .FuelCell.Sources.GasFlowBoundary                         flowCathode(
    redeclare replaceable package Medium = Medium_cath,
    X=X_air,
    m_flow=m_air,
    T=Tin_air,
    use_flow_in=true,
    use_Th_in=false)                     annotation (Placement(transformation(
          extent={{-100,-30},{-80,-10}},
                                      rotation=0)));
  .FuelCell.Sources.GasPressureBoundary                  sinkCathode(
    redeclare replaceable package Medium = Medium_cath,
    p=1.55e5,
    T=Tin_air)      annotation (Placement(transformation(extent={{100,-30},{80,
            -10}},
          rotation=0)));
  PEMFCModel.Examples.Stack coolStack(
    redeclare replaceable package Medium_an = Medium_an,redeclare replaceable package Medium_cath = Medium_cath,redeclare replaceable package Medium_cooling = MediumWater,n_cell = 11,N = 4,A_cell = 300e-4,M_stack = 10.1,L_cooling = 0.0813,D_cooling = 0.0375,p_start_in_cooling = 2.1e5,p_start_out_cooling = 1.7e5,T_start_in_cooling = 65 + 273.15,T_start_out_cooling = coolStack.T_start_in_cooling,m_flow_start_cooling = m_water,positiveFlow_cooling = true,from_dp_cooling = false,from_dp_anode = false,positiveFlow_anode = true,positiveFlow_cathode = true,from_dp_cathode = false,p_start_in_anode = 1.9e5,p_start_out_anode = 1.8e5,T_start_in_anode = Tin_fuel,T_start_out_anode = coolStack.T_start_in_anode,X_start_anode = X_fuel,m_flow_start_anode = m_fuel,p_start_in_cathode = 1.7e5,p_start_out_cathode = 1.55e5,T_start_in_cathode = Tin_air,T_start_out_cathode = coolStack.T_start_in_cathode,X_start_cathode = X_air,m_flow_start_cathode = m_air,redeclare replaceable model SubStack = .PEMFCModel.Examples.SubStack_1) annotation (Placement(transformation(extent={{-19.12918415601504,0.17416316879699423},{20.87081584398496,40.174163168796994}},rotation = 0.0,origin = {0.0,0.0})));

  .FuelCell.Sources.WaterFlowBoundary sourceW(
    redeclare replaceable package Medium = MediumWater, m_flow=8.25e-2,T = 65 + 273.15)
    annotation (Placement(transformation(extent={{-58,-50},{-38,-30}})));
  .FuelCell.Sources.WaterPressureBoundary
                         sinkP(
    redeclare replaceable package Medium = MediumWater,T = 65 + 273.15,p = 1.7e5)
    annotation (Placement(transformation(extent={{-10,10},{10,-10}},
        rotation=180,
        origin={50,-40})));
  .Modelica.Blocks.Sources.Ramp mFlow_cathode(
    duration=3600,
    startTime=0,
    height=0,
    offset=m_air)
    annotation (Placement(transformation(extent={{-109.0,-9.0},{-99.0,1.0}},rotation = 0.0,origin = {0.0,0.0})));
  .Modelica.Blocks.Sources.Ramp T_cathode(
    duration=3600,
    startTime=0,
    height=0,
    offset=Tin_air) annotation (Placement(transformation(
        extent={{5,5},{-5,-5}},
        rotation=180,
        origin={-105,11})));
  .Modelica.Blocks.Sources.Ramp mFlow_anode(
    startTime=0,
    duration=3600,
    height=0,
    offset=m_fuel)
    annotation (Placement(transformation(extent={{-110,40},{-100,50}})));
  .Modelica.Blocks.Sources.Ramp T_anode(
    duration=10,
    startTime=0,
    height=0,
    offset=Tin_fuel) annotation (Placement(transformation(
        extent={{5,5},{-5,-5}},
        rotation=180,
        origin={-105,61})));
  .Modelon.Visualizers.RealValue display_I(number=current.p.v - current.n.v,
      precision=3)
    annotation (Placement(transformation(extent={{-90,-86},{-70,-66}})));
  .Modelon.Visualizers.RealValue display_P(number=current.p.i*(current.p.v -
        current.n.v), precision=0)
    annotation (Placement(transformation(extent={{-66,-86},{-46,-66}})));
  .Modelon.Visualizers.RealValue display_P1(number=current.p.i, precision=2)
    annotation (Placement(transformation(extent={{-90,-102},{-70,-82}})));
  .Modelon.Visualizers.RealValue display_P2(number=coolStack.summary.j_external*
        1e-4*1e3, precision=1)
    annotation (Placement(transformation(extent={{-38,-86},{-18,-66}})));
  .Modelon.Visualizers.RealValue display_P3(number=sum(coolStack.coolingPipe.q.Q_flow),
      precision=1)
    annotation (Placement(transformation(extent={{-38,-102},{-18,-82}})));
  .Modelon.Visualizers.RealValue display_P7(precision=1, number=.Modelica.Units.Conversions.to_degC(coolStack.subStack.cell.T_cell_avg))
    annotation (Placement(transformation(extent={{-66,-102},{-46,-82}})));
  .FuelCell.Sensors.WaterMultiDisplaySensor
                        multiDisplaySensor4(redeclare replaceable package Medium = MediumWater)
    annotation (Placement(transformation(extent={{-36,-50},{-16,-30}})));
  .FuelCell.Sensors.MultiDisplay_phTmdot display_phTmdot4(
                                               displayUnits=true,redeclare replaceable package Medium = MediumWater)
    annotation (Placement(transformation(extent={{-36,-42},{-16,-22}})));
  .FuelCell.Sensors.MultiDisplay_phTmdot display_phTmdot1(
                                               displayUnits=true,redeclare replaceable package Medium = MediumWater)
    annotation (Placement(transformation(extent={{14,-42},{34,-22}})));
  .FuelCell.Sensors.WaterMultiDisplaySensor
                        multiDisplaySensor1(redeclare replaceable package Medium = MediumWater)
    annotation (Placement(transformation(extent={{14,-50},{34,-30}})));
  .FuelCell.Sensors.GasMultiDisplaySensor gasSensor3(redeclare replaceable package Medium = Medium_an)
    annotation (Placement(transformation(extent={{-78,20},{-58,40}})));
  .FuelCell.Sensors.MultiDisplay_phTmdot display_phTmdot3(
                                               displayUnits=true,redeclare replaceable package Medium = Medium_an)
    annotation (Placement(transformation(extent={{-78,30},{-58,50}})));
  .FuelCell.Sensors.ReformateLongCompositionDisplay display_MoleFractions1(redeclare replaceable package Medium = Medium_an)
    annotation (Placement(transformation(extent={{-54,42},{-30,62}})));
  .FuelCell.Sensors.GasMultiDisplaySensor gasSensor2(redeclare replaceable package Medium = Medium_an)
    annotation (Placement(transformation(extent={{24,20},{44,40}})));
  .FuelCell.Sensors.ReformateLongCompositionDisplay display_MoleFractions4(redeclare replaceable package Medium = Medium_an)
    annotation (Placement(transformation(extent={{50,42},{74,62}})));
  .FuelCell.Sensors.MultiDisplay_phTmdot display_phTmdot6(
                                               displayUnits=true,redeclare replaceable package Medium = Medium_an)
    annotation (Placement(transformation(extent={{24,30},{44,50}})));
  .FuelCell.Sensors.MultiDisplay_phTmdot display_phTmdot2(
                                               displayUnits=true,redeclare replaceable package Medium = Medium_cath)
    annotation (Placement(transformation(extent={{-78,-20},{-58,0}})));
  .FuelCell.Sensors.AirCompositionDisplay           display_MoleFractions3(redeclare replaceable package Medium = Medium_cath)
    annotation (Placement(transformation(extent={{-54,-14},{-32,4}})));
  .FuelCell.Sensors.GasMultiDisplaySensor gasSensor4(redeclare replaceable package Medium = Medium_cath)
    annotation (Placement(transformation(extent={{-78,-30},{-58,-10}})));
  .FuelCell.Sensors.MultiDisplay_phTmdot display_phTmdot5(
                                               displayUnits=true,redeclare replaceable package Medium = Medium_cath)
    annotation (Placement(transformation(extent={{26,-20},{46,0}})));
  .FuelCell.Sensors.GasMultiDisplaySensor gasSensor1(redeclare replaceable package Medium = Medium_cath)
    annotation (Placement(transformation(extent={{26,-30},{46,-10}})));
  .FuelCell.Sensors.AirCompositionDisplay           display_MoleFractions2(redeclare replaceable package Medium = Medium_cath)
    annotation (Placement(transformation(extent={{50,-14},{72,4}})));
initial equation
    assert(cath[1]>0,"the cathode medium does not contain the required substances",level=AssertionLevel.error);

equation
  connect(current.n,ground. p) annotation (Line(
      points={{-6,69},{-20,69},{-20,62}},
      color={0,0,255},
      smooth=Smooth.None));
  connect(current.n, coolStack.pin_n) annotation (Line(
      points={{-6,69},{-6,36.674163168796994},{-7.1291841560150395,36.674163168796994}},
      color={0,0,255},
      smooth=Smooth.None));
  connect(current.p, coolStack.pin_p) annotation (Line(
      points={{8,69},{8,36.674163168796994},{8.87081584398496,36.674163168796994}},
      color={0,0,255},
      smooth=Smooth.None));
  connect(mFlow_cathode.y, flowCathode.m_flow_in) annotation (Line(
      points={{-98.5,-4},{-96,-4},{-96,-10}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(flowCathode.T_in, T_cathode.y) annotation (Line(
      points={{-90,-10},{-90,11},{-99.5,11}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(mFlow_anode.y, flowAnode.m_flow_in) annotation (Line(
      points={{-99.5,45},{-96,45},{-96,40}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(T_anode.y, flowAnode.T_in) annotation (Line(
      points={{-99.5,61},{-90,61},{-90,40}},
      color={0,0,127},
      smooth=Smooth.None));

  connect(multiDisplaySensor4.u,display_phTmdot4. y) annotation (Line(points={{-25.95,
          -39.95},{-25.95,-34.975},{-26,-34.975},{-26,-32}}, color={0,0,0}));
  connect(multiDisplaySensor1.u,display_phTmdot1. y) annotation (Line(points={{24.05,
          -39.95},{24.05,-35.975},{24,-35.975},{24,-32}}, color={0,0,0}));
  connect(multiDisplaySensor4.portA, sourceW.fluidPort)
    annotation (Line(points={{-32,-40},{-39,-40}}, color={0,0,255}));
  connect(multiDisplaySensor4.portB, coolStack.feed_cooling)
    annotation (Line(points={{-20,-40},{-7.1291841560150395,-40},{-7.1291841560150395,3.6741631687969942}}, color={0,0,255}));
  connect(multiDisplaySensor1.portB, sinkP.fluidPort)
    annotation (Line(points={{30,-40},{41,-40}}, color={0,0,255}));
  connect(multiDisplaySensor1.portA, coolStack.drain_cooling)
    annotation (Line(points={{18,-40},{8.87081584398496,-40},{8.87081584398496,3.6741631687969942}}, color={0,0,255}));
  connect(display_phTmdot3.y,gasSensor3. u) annotation (Line(points={{-68,40},{
          -68,30.05},{-67.95,30.05}},
                                  color={0,0,0}));
  connect(display_MoleFractions1.data,gasSensor3. u) annotation (Line(points={{-42,37.2941},{-42,32},{-67.95,32},{-67.95,30.05}},
                                                         color={0,0,0}));
  connect(gasSensor3.portA, flowAnode.fluidPort)
    annotation (Line(points={{-74,30},{-81,30}}, color={255,128,0}));
  connect(gasSensor3.portB, coolStack.feed_an)
    annotation (Line(points={{-62,30},{-17.12918415601504,30},{-17.12918415601504,30.174163168796994}}, color={255,128,0}));
  connect(gasSensor2.u,display_phTmdot6. y) annotation (Line(points={{34.05,
          30.05},{34.05,35.025},{34,35.025},{34,40}},
                                               color={0,0,0}));
  connect(gasSensor2.u,display_MoleFractions4. data) annotation (Line(points={{34.05,30.05},{34.05,32},{62,32},{62,37.2941}},
                                                   color={0,0,0}));
  connect(gasSensor2.portA, coolStack.drain_an)
    annotation (Line(points={{28,30},{18.87081584398496,30},{18.87081584398496,30.174163168796994}}, color={255,128,0}));
  connect(gasSensor2.portB, sinkAnode.fluidPort)
    annotation (Line(points={{40,30},{81,30}}, color={255,128,0}));
  connect(gasSensor4.u,display_phTmdot2. y) annotation (Line(points={{-67.95,
          -19.95},{-67.95,-14.975},{-68,-14.975},{-68,-10}},
                                                     color={0,0,0}));
  connect(gasSensor4.u,display_MoleFractions3. data) annotation (Line(points={{-67.95,
          -19.95},{-67.95,-17.975},{-43,-17.975},{-43,-14}}, color={0,0,0}));
  connect(gasSensor4.portA, flowCathode.fluidPort)
    annotation (Line(points={{-74,-20},{-81,-20}}, color={255,128,0}));
  connect(gasSensor4.portB, coolStack.feed_cath) annotation (Line(points={{-62,-20},{-17.12918415601504,-20},{-17.12918415601504,10.174163168796994}}, color={255,128,0}));
  connect(gasSensor1.u,display_MoleFractions2. data) annotation (Line(points={{36.05,
          -19.95},{36.05,-17.975},{61,-17.975},{61,-14}},    color={0,0,0}));
  connect(sinkCathode.fluidPort, gasSensor1.portB)
    annotation (Line(points={{81,-20},{42,-20}}, color={255,128,0}));
  connect(gasSensor1.portA, coolStack.drain_cath)
    annotation (Line(points={{30,-20},{18.87081584398496,-20},{18.87081584398496,10.174163168796994}}, color={255,128,0}));
  connect(gasSensor1.u, display_phTmdot5.y) annotation (Line(points={{36.05,
          -19.95},{36.05,-14.975},{36,-14.975},{36,-10}}, color={0,0,0}));
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
end CoolStack_original;
