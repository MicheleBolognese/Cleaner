within PEMFCModel.Stacks.Interfaces;
partial model LumpedCathodeGeometry "Cathode channel geometry parameters"

  constant Real pi = .Modelica.Constants.pi;
  parameter Boolean enable_setting_cathode_geometry = true "If true, selectability of cathode parameters enabled" annotation(Evaluate=true,Dialog(tab = "Substack",group="Cathode side geometry"));
  
  parameter Real n_channels_cathode(min = 1.0) = 1.0 "Number of parallel channels at cathode" annotation(Evaluate=true,Dialog(tab = "Substack", group="Cathode side geometry", enable = enable_setting_cathode_geometry));
  parameter .Modelica.Units.SI.Length L_cathode = 0.1 "Length of cathode channel"  annotation(Evaluate=true,Dialog(tab = "Substack", group="Cathode side geometry", enable = enable_setting_cathode_geometry));
  parameter .Modelica.Units.SI.Length D_cathode = 0.01 "Diameter of single cathode channel" annotation(Evaluate=true,Dialog(tab = "Substack", group="Cathode side geometry", enable = enable_setting_cathode_geometry));
  parameter .Modelica.Units.SI.Length C_cathode = pi*D_cathode "Circumference of single cathode channel" annotation(Evaluate=true,Dialog(tab = "Substack", group="Cathode side geometry", enable = enable_setting_cathode_geometry));
  parameter .Modelica.Units.SI.Area A_cathode = pi*D_cathode*D_cathode/4 "Cross section area of single cathode channel" annotation(Evaluate=true,Dialog(tab = "Substack", group="Cathode side geometry", enable = enable_setting_cathode_geometry));
  parameter .Modelica.Units.SI.Volume V_cathode = n_channels_cathode*L_cathode*A_cathode "Total volume of cathode channels" annotation(Evaluate=true,Dialog(tab = "Substack", group="Cathode side geometry", enable = enable_setting_cathode_geometry));
  parameter .Modelica.Units.SI.Area A_heat_cathode = C_cathode*L_cathode "Heat transfer area of single cathode channel" annotation(Evaluate=true,Dialog(tab = "Substack", group="Cathode side geometry", enable = enable_setting_cathode_geometry));

  final parameter .Modelica.Units.SI.Length Dhyd_cathode = .Modelon.ThermoFluid.Functions.CharacteristicNumbers.HydraulicDiameter(A_cathode, C_cathode) "Hydraulic diameter of single cathode channel" annotation (Dialog(group="Cathode side geometry",enable = true));

  annotation (Documentation(revisions="<html>
Copyright &copy; 2004-2024, MODELON AB <br /> The use of this software component is regulated by the licensing conditions for Modelon Libraries. <br />This copyright notice must, unaltered, accompany all components that are derived from, copied from, <br />or by other means have their origin from any Modelon Library.
</html>", info="<html>
<h4>LumpedCathodeGeometry</h4>
<p>Interface with cathode channel geometry parameters.</p>
</html>"));

end LumpedCathodeGeometry;