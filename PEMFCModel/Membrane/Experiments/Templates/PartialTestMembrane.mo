within PEMFCModel.Membrane.Experiments.Templates;
partial model PartialTestMembrane "Test equivalent circuit cell"

  extends .Modelon.Icons.Experiment;

  replaceable package Medium_an = .FuelCell.Media.PreDefined.IdealGases.NASAReformateLong constrainedby
    .FuelCell.Media.Templates.ReactionGas;
  replaceable package Medium_cath = .FuelCell.Media.PreDefined.IdealGases.NASAMoistAir constrainedby
    .FuelCell.Media.Templates.ReactionGas;

  // Experiment conditions
  parameter .Modelica.Units.SI.Pressure p_anode = 1.013e5;
  parameter .Modelica.Units.SI.Pressure p_cathode = p_anode;
  parameter .Modelica.Units.SI.Temperature T_anode = 353.15;
  parameter .Modelica.Units.SI.Temperature T_cathode = T_anode;
  parameter .Modelica.Units.SI.MassFraction X_anode[Medium_an.nS] = Medium_an.moleToMassFractions({0.33, 0.1, 0.06, 0.06, 0.2, 0.25, 0}, Medium_an.MMX);
  parameter .Modelica.Units.SI.MassFraction X_cathode[Medium_cath.nS] = Medium_cath.reference_X;
  
  parameter .Modelica.Units.SI.SpecificEnthalpy h_anode = Medium_an.specificEnthalpy_pTX(p_anode, T_anode, X_anode);

  parameter .Modelica.Units.SI.SpecificEnthalpy h_cathode = Medium_cath.specificEnthalpy_pTX(p_anode, T_anode, X_cathode);

  replaceable PEMFCModel.Membrane.Templates.PartialCell cellMembrane(
    redeclare package Medium_an = Medium_an,
    redeclare package Medium_cath = Medium_cath,
    C_cell = 250,
    N = 1,
    Tstart = T_anode,
    pstart = p_anode,
    A_cell = 361e-4,
    n_cell = 50)
  constrainedby PEMFCModel.Membrane.Templates.PartialCell(
    redeclare package Medium_an = Medium_an,
    redeclare package Medium_cath = Medium_cath,
    C_cell = 250,
    N = 1,
    Tstart = T_anode,
    pstart = p_anode,
    A_cell = 361e-4,
    n_cell = 50) annotation(Placement(transformation(extent={{-20,-20},{20,20}}, rotation=0)));
  
  .Modelica.Thermal.HeatTransfer.Sources.PrescribedTemperature temp_cath[cellMembrane.N] annotation(Placement(transformation(
        origin={40,-30},
        extent={{-10,-10},{10,10}},
        rotation=180)));
  
  .Modelica.Thermal.HeatTransfer.Sources.PrescribedTemperature temp_anode[cellMembrane.N] annotation(Placement(transformation(
        origin={40,30},
        extent={{-10,-10},{10,10}},
        rotation=180)));
  
  .Modelica.Blocks.Sources.Ramp ramp_cath[cellMembrane.N](
    each startTime = 200,
    each duration = 300,
    each height = 0,
    each offset = T_cathode) annotation(Placement(transformation(extent={{90,-40},{70,-20}}, rotation=0)));

   .Modelica.Blocks.Sources.Ramp ramp_anode[cellMembrane.N](
    each startTime = 200,
    each duration = 300,
    each height = 0,
    each offset = T_anode) annotation(Placement(transformation(extent={{90,20},{70,40}},  rotation=0)));

  .FuelCell.Sources.FixedMassBoundary anode[cellMembrane.N](
    redeclare package Medium = Medium_an,
    each p = p_anode,
    each h = h_anode,
    each X = X_anode) annotation(Placement(transformation(extent={{-102.0,30.0},{-82.0,50.0}}, rotation=0.0,origin = {0.0,0.0})));
  
  .FuelCell.Sources.FixedMassBoundary cathode[cellMembrane.N](
    redeclare package Medium = Medium_cath,
    each p = p_cathode,
    each h = h_cathode,
    each X = X_cathode) annotation(Placement(transformation(extent={{-102,-50},{-82,-30}}, rotation=0)));
  
  .Modelica.Electrical.Analog.Basic.Ground ground annotation(Placement(transformation(extent={{-38.14853493147929,8.899322936014801},{-26.14853493147929,20.8993229360148}},
                                                                                                                   rotation=0.0,origin = {0.0,0.0})));
  
  .Modelica.Electrical.Analog.Sources.RampCurrent current(
    duration = 200,
    offset = 0,
    I = 150) annotation(Placement(transformation(
        origin={-40.98143313356508,0.3271443778550225},
        extent={{-10.0,10.0},{10.0,-10.0}},
        rotation=90.0)));
  
  .Modelon.Visualizers.RealValue display_V(number = current.p.v - current.n.v, precision = 1) annotation(Placement(transformation(extent={{-38,70},{-18,90}})));
  
  .Modelon.Visualizers.RealValue display_P(number = current.p.i*(current.p.v - current.n.v), precision = 0) annotation(Placement(transformation(extent={{-14,70},{6,90}})));

  .Modelon.Visualizers.RealValue display_T(number = .Modelica.Units.Conversions.to_degC(cellMembrane.T_cell_avg), precision = 1) annotation(Placement(transformation(extent={{10,70},{30,90}})));

  .FuelCell.Sensors.GasMultiDisplaySensor_MassPort gasMultiDisplaySensor_an[cellMembrane.N](redeclare package Medium = Medium_an) annotation(Placement(transformation(extent={{-72,50},{-52,30}})));

  .FuelCell.Sensors.GasMultiDisplaySensor_MassPort gasMultiDisplaySensor_cath[cellMembrane.N](redeclare package Medium = Medium_cath) annotation(Placement(transformation(extent={{-72,-50},{-52,-30}})));

  .FuelCell.Sensors.AirCompositionDisplay airCompositionDisplay_cath(
    redeclare package Medium = Medium_cath,
    sensorType = .FuelCell.Sensors.Types.SensorType.Flow,
    precision = 4,
    displayMassUnit = true) annotation (Placement(transformation(extent={{-76,-34},{-48,-10}})));
  
  .FuelCell.Sensors.MultiDisplay_phTmdot multiDisplay_phTmdot_cath(redeclare package Medium = Medium_cath, displayUnits = true) annotation(Placement(transformation(extent={{-74,-64},{-50,-40}})));
  
  .FuelCell.Sensors.ReformateLongCompositionDisplay reformateLongCompositionDisplay_an(
    redeclare package Medium = Medium_an,
    sensorType = .FuelCell.Sensors.Types.SensorType.Flow,
    precision = 4,
    displayMassUnit = true) annotation(Placement(transformation(extent={{-76,28},{-48,4}})));
  
  .FuelCell.Sensors.MultiDisplay_phTmdot multiDisplay_phTmdot_an(displayUnits = true, redeclare package Medium = Medium_an) annotation(Placement(transformation(extent={{-74,40},{-50,64}})));

equation

  connect(cellMembrane.pin_p, current.p) annotation(Line(
      points={{-18,-10},{-40.98143313356508,-10},{-40.98143313356508,-9.672855622144978}},
      color={0,0,255},
      smooth=Smooth.None));
  connect(cellMembrane.pin_n, current.n) annotation(Line(
      points={{-17.937388637815914,10},{-40.98143313356508,10},{-40.98143313356508,10.327144377855022}},
      color={0,0,255},
      smooth=Smooth.None));
  connect(current.n, ground.p) annotation(Line(
      points={{-40.98143313356508,10.327144377855022},{-40.98143313356508,20.8993229360148},{-32.14853493147929,20.8993229360148}},
      color={0,0,255},
      smooth=Smooth.None));

  connect(ramp_anode.y, temp_anode.T) annotation(Line(points={{69,30},{52,30}}, color={0,0,127}));
  connect(ramp_cath.y, temp_cath.T) annotation(Line(points={{69,-30},{52,-30}}, color={0,0,127}));
  connect(temp_cath.port, cellMembrane.wall_cath) annotation(Line(
      points={{30,-30},{8,-30},{8,-10}},
      color={191,0,0},
      smooth=Smooth.None));
  connect(cellMembrane.wall_an, temp_anode.port) annotation(Line(points={{8,10},{8,30},{30,30}}, color={191,0,0}));
  
  connect(gasMultiDisplaySensor_an.portB, cellMembrane.port_an) annotation(Line(points={{-56,40},{-6,40},{-6,8.8}},color={0,0,255}));
  connect(gasMultiDisplaySensor_cath.portB, cellMembrane.port_cath) annotation(Line(points={{-56,-40},{-6,-40},{-6,-8.8}},
                                                          color={0,0,255}));
  connect(airCompositionDisplay_cath.data, gasMultiDisplaySensor_cath[1].u) annotation(Line(points={{-62,-34},{-62,-37.975},{-62.05,-37.975},{-62.05,-39.95}},
                                                                      color={0,0,0}));
  connect(cathode.port, gasMultiDisplaySensor_cath.portA) annotation(Line(points={{-82.6,-40},{-68,-40}},           color={0,0,255}));
  connect(anode.port, gasMultiDisplaySensor_an.portA) annotation(Line(points={{-82.6,40},{-68,40}},          color={0,0,255}));
  connect(reformateLongCompositionDisplay_an.data, gasMultiDisplaySensor_an[1].u) annotation(Line(points={{-62,33.6471},{-62,34.7985},{-62.05,34.7985},{-62.05,39.95}},
                                                                        color={0,0,0})); 
  connect(multiDisplay_phTmdot_an.y, gasMultiDisplaySensor_an[1].u) annotation(Line(points={{-62,52},{-62,46.975},{-62.05,46.975},{-62.05,39.95}}, color={0,0,0}));
  connect(multiDisplay_phTmdot_cath.y, gasMultiDisplaySensor_cath[1].u) annotation(Line(points={{-62,-52},{-62,-44.975},{-62.05,-44.975},{-62.05,-39.95}}, color={0,0,0}));
  
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
Copyright &copy; 2004-2024, MODELON AB <br /> The use of this software component is regulated by the licensing conditions for Modelon Libraries. <br />This copyright notice must, unaltered, accompany all components that are derived from, copied from, <br />or by other means have their origin from any Modelon Library.
</html>", info="<html>
<h4>PartialTestMembrane</h4>
<p>Experiment setup for testing cell membranes with a linearly increasing current.</p>
</html>"));

end PartialTestMembrane;