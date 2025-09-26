within SAMU;

partial model PartialTestumi "Test equivalent circuit cell"

  extends .Modelon.Icons.Experiment;
  parameter .Modelica.Units.SI.Pressure p_anode=1.013e5;
  parameter .Modelica.Units.SI.Pressure p_cathode=p_anode;
  parameter .Modelica.Units.SI.Temperature T_anode=1033;
  parameter .Modelica.Units.SI.Temperature T_cathode=T_anode;
  parameter .Modelica.Units.SI.MassFraction[Medium_an.nS] X_anode=Medium_an.moleToMassFractions({0.33,0.1,0.06,0.06,0.2,0.25,
      0}, Medium_an.MMX);
  parameter .Modelica.Units.SI.MassFraction[Medium_cath.nS] X_cathode=Medium_cath.reference_X[1:Medium_cath.nS];

  // Experiment conditions
  replaceable package Medium_an =
      .FuelCell.Media.PreDefined.IdealGases.NASAReformateLong                             constrainedby
    .FuelCell.Media.Templates.ReactionGas;
  replaceable package Medium_cath =
      .FuelCell.Media.PreDefined.IdealGases.NASAMoistAir                               constrainedby
    .FuelCell.Media.Templates.ReactionGas;

  parameter .Modelica.Units.SI.SpecificEnthalpy h_anode=Medium_an.specificEnthalpy_pTX(
      p_anode,
      T_anode,
      anode[1].X);

  parameter .Modelica.Units.SI.SpecificEnthalpy h_cathode=Medium_cath.specificEnthalpy_pTX(
      p_anode,
      T_anode,
      cathode[1].X);
  .Modelica.Thermal.HeatTransfer.Sources.PrescribedTemperature temp_cath[cellMembrane.N]
    annotation (Placement(transformation(
        origin={40,-30},
        extent={{-10,-10},{10,10}},
        rotation=180)));

  .FuelCell.Sources.FixedMassBoundary anode[cellMembrane.N](
    redeclare package Medium = Medium_an,
    each p=p_anode,
    each h=h_anode,
    each X=X_anode) annotation (Placement(transformation(extent={{-102,30},{-82,50}}, rotation=0)));
  .FuelCell.Sources.FixedMassBoundary cathode[cellMembrane.N](
    redeclare package Medium = Medium_cath,
    each p=p_cathode,
    each h=h_cathode,
    each X=X_cathode) annotation (Placement(transformation(extent={{-102,-50},{-82,-30}}, rotation=0)));
  .Modelica.Blocks.Sources.Ramp ramp_anode[cellMembrane.N](
    each startTime=200,
    each duration=300,
    each height=0,
    each offset=T_anode) annotation (Placement(transformation(extent={{90,20},{70,40}},  rotation=0)));
  .Modelon.Visualizers.RealValue display_I( precision=1)
    annotation (Placement(transformation(extent={{-38,70},{-18,90}})));
  .Modelon.Visualizers.RealValue display_P( precision=0)
    annotation (Placement(transformation(extent={{-14.0,70.0},{6.0,90.0}},rotation = 0.0,origin = {0.0,0.0})));
  .Modelon.Visualizers.RealValue display_T(
                                                                      precision=1)
    annotation (Placement(transformation(extent={{10.0,70.0},{30.0,90.0}},rotation = 0.0,origin = {0.0,0.0})));
  .FuelCell.Sensors.GasMultiDisplaySensor_MassPort gasMultiDisplaySensor_an[cellMembrane.N](redeclare package Medium = Medium_an)
    annotation (Placement(transformation(extent={{-72,50},{-52,30}})));
  .FuelCell.Sensors.GasMultiDisplaySensor_MassPort gasMultiDisplaySensor_cath[cellMembrane.N](redeclare package Medium =
        Medium_cath)
    annotation (Placement(transformation(extent={{-72,-50},{-52,-30}})));
  .FuelCell.Sensors.AirCompositionDisplay airCompositionDisplay_cath(
    redeclare package Medium = Medium_cath,
    sensorType=.FuelCell.Sensors.Types.SensorType.Flow,
    precision=4,
    displayMassUnit=true) annotation (Placement(transformation(extent={{-76,-34},{-48,-10}})));
  .FuelCell.Sensors.MultiDisplay_phTmdot multiDisplay_phTmdot_cath(redeclare package Medium =
        Medium_cath,                                                                             displayUnits=true)
    annotation (Placement(transformation(extent={{-74,-64},{-50,-40}})));
  .FuelCell.Sensors.ReformateLongCompositionDisplay reformateLongCompositionDisplay_an(
    redeclare package Medium = Medium_an,
    sensorType=.FuelCell.Sensors.Types.SensorType.Flow,
    precision=4,
    displayMassUnit=true) annotation (Placement(transformation(extent={{-76,28},{-48,4}})));
  .FuelCell.Sensors.MultiDisplay_phTmdot multiDisplay_phTmdot_an(displayUnits=true, redeclare package Medium = Medium_an)
    annotation (Placement(transformation(extent={{-74,40},{-50,64}})));

  .Modelica.Thermal.HeatTransfer.Sources.PrescribedTemperature temp_anode[cellMembrane.N]
    annotation (Placement(transformation(
        origin={40,30},
        extent={{-10,-10},{10,10}},
        rotation=180)));
  .Modelica.Blocks.Sources.Ramp ramp_cath[cellMembrane.N](
    each startTime=200,
    each duration=300,
    each height=0,
    each offset=T_cathode) annotation (Placement(transformation(extent={{90,-40},{70,-20}}, rotation=0)));
    replaceable SAMU.PartialCellTransport_humidifier cellMembrane constrainedby SAMU.PartialCellTransport_humidifier annotation(Placement(transformation(extent = {{-24,-22},{16,18}},rotation = 0,origin = {0,0})));

equation
  connect(airCompositionDisplay_cath.data, gasMultiDisplaySensor_cath[1].u)
    annotation (Line(points={{-62,-34},{-62,-37.975},{-62.05,-37.975},{-62.05,-39.95}},
                                                                      color={0,0,0}));
  connect(cathode.port, gasMultiDisplaySensor_cath.portA)
    annotation (Line(points={{-82.6,-40},{-68,-40}},           color={0,0,255}));
  connect(anode.port, gasMultiDisplaySensor_an.portA) annotation (Line(points={{-82.6,40},{-68,40}},          color={0,0,255}));
  connect(reformateLongCompositionDisplay_an.data, gasMultiDisplaySensor_an[1].u)
    annotation (Line(points={{-62,33.6471},{-62,34.7985},{-62.05,34.7985},{-62.05,39.95}},
                                                                        color={0,0,0}));
  connect(ramp_anode.y, temp_anode.T) annotation (Line(points={{69,30},{52,30}}, color={0,0,127}));
  connect(ramp_cath.y, temp_cath.T) annotation (Line(points={{69,-30},{52,-30}}, color={0,0,127}));
  connect(multiDisplay_phTmdot_an.y, gasMultiDisplaySensor_an[1].u)
    annotation (Line(points={{-62,52},{-62,46.975},{-62.05,46.975},{-62.05,39.95}}, color={0,0,0}));
  connect(multiDisplay_phTmdot_cath.y, gasMultiDisplaySensor_cath[1].u)
    annotation (Line(points={{-62,-52},{-62,-44.975},{-62.05,-44.975},{-62.05,-39.95}}, color={0,0,0}));
    connect(temp_cath.port,cellMembrane.wall_cath) annotation(Line(points = {{30,-30},{4,-30},{4,-12}},color = {191,0,0}));
    connect(gasMultiDisplaySensor_an.portB,cellMembrane.port_an) annotation(Line(points = {{-56,40},{-10,40},{-10,6.8}},color = {0,0,255}));
    connect(gasMultiDisplaySensor_cath.portB,cellMembrane.port_cath) annotation(Line(points = {{-56,-40},{-10,-40},{-10,-10.8}},color = {0,0,255}));
    connect(cellMembrane.wall_an,temp_anode.port) annotation(Line(points = {{4,8},{4,30},{30,30}},color = {191,0,0}));
  annotation (
    Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{100,100}}), graphics={
        Rectangle(
          extent={{-40,98},{40,74}},
          lineColor={215,215,215},
          fillColor={215,215,215},
          fillPattern=FillPattern.Solid,
          radius=2),
        Text(
          extent={{-36,90},{-18,84}},
          lineColor={0,0,0},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid,
          textString="Voltage [V]"),
        Text(
          extent={{-14,90},{4,84}},
          lineColor={0,0,0},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid,
          textString="Power [W]"),
        Line(points={{-40,90},{40,90}}, color={0,0,0}),
        Text(
          extent={{-10,90},{8,98}},
          lineColor={0,0,0},
          textStyle={TextStyle.Bold},
          textString="Cell"),
        Text(
          extent={{6,90},{38,84}},
          lineColor={0,0,0},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid,
          textString="Temperature [degC]")}),
    Documentation(revisions="<html>
Copyright &copy; 2004-2025, MODELON AB <br /> The use of this software component is regulated by the licensing conditions for Modelon Libraries. <br />This copyright notice must, unaltered, accompany all components that are derived from, copied from, <br />or by other means have their origin from any Modelon Library.
</html>", info="<html>
<h4>PartialTestMembrane</h4>
<p>Experiment setup for testing cell membranes with a linearly increasing current.</p>
</html>"));
end PartialTestumi;
