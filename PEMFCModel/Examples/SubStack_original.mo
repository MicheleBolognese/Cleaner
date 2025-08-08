within PEMFCModel.Examples;

model SubStack_original "PEMFC substack test model with prox at anode inlet"
  extends .Modelon.Icons.Experiment;

  replaceable package Medium_an =
      .FuelCell.Media.PreDefined.IdealGases.NASAReformateShort                             constrainedby
    .FuelCell.Media.Templates.ReactionGas "reformate medium";
  replaceable package Medium_cath =
      .FuelCell.Media.PreDefined.IdealGases.NASAMoistAir                               constrainedby
    .FuelCell.Media.Templates.ReactionGas "air medium";

  .FuelCell.Sources.GasFlowBoundary                      flowCathode(
    redeclare package Medium = Medium_cath,
    m_flow=4.711e-2,
    T=333.15)                          annotation (Placement(transformation(
          extent={{-100,-20},{-80,0}},rotation=0)));

  .FuelCell.Stacks.PEMFC.SubStack subStack(
    redeclare package Medium_an = Medium_an,
    redeclare package Medium_cath = Medium_cath,redeclare replaceable model Membrane = .FuelCell.Membranes.PEMFC.Empirical,
    N=4,
    n_cell=1,
    redeclare model Friction_anode =
        .Modelon.ThermoFluid.FlowChannels.PipeResistances.SinglePhase.LinearOperatingPointLoss
        (
        d0=0.6,
        dp0(displayUnit="kPa") = 500,
        m_flow0=2e-6),
    redeclare model Friction_cathode =
        .Modelon.ThermoFluid.FlowChannels.PipeResistances.SinglePhase.LinearOperatingPointLoss
        (
        d0=0.5,
        dp0(displayUnit="kPa") = 500,
        m_flow0=9e-5),
    redeclare model HeatTransfer_cathode =
        .Modelon.ThermoFluid.FlowChannels.HeatTransfer.SinglePhase.ConstantCoefficient
        (                                                                                                                 alpha=
            100),
    redeclare model HeatTransfer_anode =
        .Modelon.ThermoFluid.FlowChannels.HeatTransfer.SinglePhase.ConstantCoefficient
        (                                                                                                               alpha=100),addProxToAnode = false,M_cell = 30 / 50)
    annotation (Placement(transformation(extent={{-20,0},{20,40}}, rotation=0)));

  .FuelCell.Sources.GasFlowBoundary                      flowAnode(
    redeclare package Medium = Medium_an,
    m_flow=1.035e-3,
    X=Medium_an.moleToMassFractions({0.48,0.006,0.07,0.162,0.162,0.12},
        Medium_an.MMX),
    T=333.15)                          annotation (Placement(transformation(
          extent={{-100,18},{-80,38}},
                                     rotation=0)));

  .FuelCell.Sources.GasPressureBoundary             sinkCathode(
    redeclare package Medium = Medium_cath,
    p=101300,
    T=1073.15)      annotation (Placement(transformation(extent={{98,-20},{78,0}},
          rotation=0)));
  .FuelCell.Sources.GasPressureBoundary             sinkAnode(
    redeclare package Medium = Medium_an,
    p=101300,
    T=1073.15)      annotation (Placement(transformation(extent={{98,18},{78,38}},
          rotation=0)));

  .Modelica.Electrical.Analog.Basic.Ground ground
    annotation (Placement(transformation(extent={{-16,-52},{0,-36}},   rotation=
           0)));
  .Modelica.Electrical.Analog.Sources.RampCurrent current(
    startTime=0,
    offset=10,
    I=320,
    duration=4000)
    annotation (Placement(transformation(
        origin={-18,-26},
        extent={{-10,-10},{10,10}},
        rotation=0)));
  .Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow cooling[subStack.N](each alpha=0.25, each T_ref=
        .Modelica.Units.Conversions.from_degC(800)) annotation (Placement(transformation(
        origin={8,52},
        extent={{-10,-10},{10,10}},
        rotation=270)));
  .Modelica.Blocks.Sources.Ramp Q_flow[subStack.N](
    each startTime=20,
    each duration=50,
    each height=-0/subStack.N) annotation (Placement(transformation(extent={{-18,64},
            {-2,80}},         rotation=0)));
  .Modelon.Visualizers.RealValue display_I(number=current.p.v - current.n.v,
      precision=3)
    annotation (Placement(transformation(extent={{-90,-86},{-70,-66}})));
  .Modelon.Visualizers.RealValue display_P(number=current.p.i*(current.p.v -
        current.n.v), precision=0)
    annotation (Placement(transformation(extent={{-66,-86},{-46,-66}})));
  .Modelon.Visualizers.RealValue display_P1(number=current.p.i, precision=2)
    annotation (Placement(transformation(extent={{-90,-102},{-70,-82}})));
  .Modelon.Visualizers.RealValue display_P2(precision=1, number=subStack.summary.j_external
        *1e-4*1e3)
    annotation (Placement(transformation(extent={{-38,-86},{-18,-66}})));
  .Modelon.Visualizers.RealValue display_P3(precision=1, number=-sum(subStack.wall.Q_flow))
    annotation (Placement(transformation(extent={{-38,-102},{-18,-82}})));
  .Modelon.Visualizers.RealValue display_P7(precision=1, number=.Modelica.Units.Conversions.to_degC(subStack.cell.T_cell_avg))
    annotation (Placement(transformation(extent={{-66,-102},{-46,-82}})));
  .FuelCell.Sensors.MultiDisplay_phTmdot display_phTmdot3(
                                               displayUnits=true, redeclare
      package                                                                       Medium =
                       Medium_an)
    annotation (Placement(transformation(extent={{-76,28},{-56,48}})));
  .FuelCell.Sensors.ReformateLongCompositionDisplay display_MoleFractions1(redeclare
      package                                                                      Medium =
                       Medium_an)
    annotation (Placement(transformation(extent={{-52,40},{-28,60}})));
  .FuelCell.Sensors.GasMultiDisplaySensor gasSensor3(redeclare package Medium = Medium_an)
    annotation (Placement(transformation(extent={{-76,18},{-56,38}})));
  .FuelCell.Sensors.MultiDisplay_phTmdot display_phTmdot2(
                                               displayUnits=true, redeclare
      package                                                                       Medium =
                       Medium_cath)
    annotation (Placement(transformation(extent={{-76,-10},{-56,10}})));
  .FuelCell.Sensors.GasMultiDisplaySensor gasSensor4(redeclare package Medium =
        Medium_cath)
    annotation (Placement(transformation(extent={{-76,-20},{-56,0}})));
  .FuelCell.Sensors.AirCompositionDisplay           display_MoleFractions3(redeclare
      package                                                                      Medium =
                       Medium_cath)
    annotation (Placement(transformation(extent={{-52,-4},{-30,14}})));
  .FuelCell.Sensors.MultiDisplay_phTmdot display_phTmdot6(
                                               displayUnits=true, redeclare
      package                                                                       Medium =
                       Medium_an)
    annotation (Placement(transformation(extent={{28,26},{48,46}})));
  .FuelCell.Sensors.GasMultiDisplaySensor gasSensor2(redeclare package Medium = Medium_an)
    annotation (Placement(transformation(extent={{28,18},{48,38}})));
  .FuelCell.Sensors.ReformateLongCompositionDisplay display_MoleFractions4(redeclare
      package                                                                      Medium =
                       Medium_an)
    annotation (Placement(transformation(extent={{52,38},{76,58}})));
  .FuelCell.Sensors.MultiDisplay_phTmdot display_phTmdot5(
                                               displayUnits=true, redeclare
      package                                                                       Medium =
                       Medium_cath)
    annotation (Placement(transformation(extent={{32,-12},{52,8}})));
  .FuelCell.Sensors.GasMultiDisplaySensor gasSensor1(redeclare package Medium =
        Medium_cath)
    annotation (Placement(transformation(extent={{32,-20},{52,0}})));
  .FuelCell.Sensors.AirCompositionDisplay           display_MoleFractions2(redeclare
      package                                                                      Medium =
                       Medium_cath)
    annotation (Placement(transformation(extent={{56,-6},{76,12}})));
equation
  connect(current.n,ground. p) annotation (Line(
      points={{-8,-26},{-8,-36}},
      color={0,0,255},
      smooth=Smooth.None));
  connect(current.n, subStack.pin_n) annotation (Line(
      points={{-8,-26},{0,-26},{0,6.8}},
      color={0,0,255},
      smooth=Smooth.None));
  connect(Q_flow.y, cooling.Q_flow) annotation (Line(
      points={{-1.2,72},{8.25,72},{8.25,62},{8,62}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(cooling.port, subStack.wall) annotation (Line(
      points={{8,42},{8,20}},
      color={191,0,0},
      smooth=Smooth.None));

  connect(current.p, subStack.pin_p) annotation (Line(points={{-28,-26},{-28,40},
          {0,40},{0,33.2}},               color={0,0,255}));
  connect(display_phTmdot3.y,gasSensor3. u) annotation (Line(points={{-66,38},{
          -66,28.05},{-65.95,28.05}},
                                  color={0,0,0}));
  connect(display_MoleFractions1.data,gasSensor3. u) annotation (Line(points={{-40,
          35.2941},{-40,30},{-65.95,30},{-65.95,28.05}}, color={0,0,0}));
  connect(gasSensor3.portA, flowAnode.fluidPort)
    annotation (Line(points={{-72,28},{-81,28}}, color={255,128,0}));
  connect(gasSensor3.portB, subStack.feed_anode)
    annotation (Line(points={{-60,28},{-18,28}}, color={255,128,0}));
  connect(gasSensor4.u,display_phTmdot2. y) annotation (Line(points={{-65.95,
          -9.95},{-65.95,-4.975},{-66,-4.975},{-66,0}},
                                                     color={0,0,0}));
  connect(gasSensor4.u,display_MoleFractions3. data) annotation (Line(points={{-65.95,
          -9.95},{-65.95,-7.975},{-41,-7.975},{-41,-4}},     color={0,0,0}));
  connect(subStack.feed_cathode, gasSensor4.portB) annotation (Line(points={{
          -18,12},{-18,-10},{-60,-10}}, color={255,128,0}));
  connect(gasSensor4.portA, flowCathode.fluidPort)
    annotation (Line(points={{-72,-10},{-81,-10}}, color={255,128,0}));
  connect(display_MoleFractions4.data, gasSensor2.u) annotation (Line(points={{64,
          33.2941},{64,30},{38,30},{38,28.05},{38.05,28.05}},    color={0,0,0}));
  connect(gasSensor2.u, display_phTmdot6.y) annotation (Line(points={{38.05,
          28.05},{38.05,32.025},{38,32.025},{38,36}}, color={0,0,0}));
  connect(sinkAnode.fluidPort, gasSensor2.portB)
    annotation (Line(points={{79,28},{44,28}}, color={255,128,0}));
  connect(gasSensor2.portA, subStack.drain_anode)
    annotation (Line(points={{32,28},{18,28}}, color={255,128,0}));
  connect(gasSensor1.u,display_phTmdot5. y) annotation (Line(points={{42.05,
          -9.95},{42.05,-4.975},{42,-4.975},{42,-2}},
                                                  color={0,0,0}));
  connect(gasSensor1.u, display_MoleFractions2.data) annotation (Line(points={{
          42.05,-9.95},{42.05,-8},{66,-8},{66,-6}}, color={0,0,0}));
  connect(sinkCathode.fluidPort, gasSensor1.portB)
    annotation (Line(points={{79,-10},{48,-10}}, color={255,128,0}));
  connect(gasSensor1.portA, subStack.drain_cathode)
    annotation (Line(points={{36,-10},{18,-10},{18,12}}, color={255,128,0}));
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
end SubStack_original;
