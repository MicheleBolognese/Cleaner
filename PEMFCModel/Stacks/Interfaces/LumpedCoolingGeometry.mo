within PEMFCModel.Stacks.Interfaces;
partial model LumpedCoolingGeometry "Geometry parameters for lumped cooling channel or pipe"
  
  constant Real pi = .Modelica.Constants.pi;
  
  parameter Boolean enable_setting_cooling_geometry = true "If true, selectability of cooling channel initialization parameter enabled" annotation(Dialog(tab = "Cooling",group="Geometry"));
  parameter Real n_channels_cooling(min = 1.0) = 1.0 "Number of parallel channels" annotation(Dialog(tab = "Cooling",group="Geometry",enable = enable_setting_cooling_geometry ));
  parameter .Modelica.Units.SI.Length L_cooling = 0.1 "Cooling channel length" annotation (Dialog(tab = "Cooling",group="Geometry",enable = enable_setting_cooling_geometry));
  parameter .Modelica.Units.SI.Length D_cooling = 0.01 "Diameter of single cooling channel" annotation (Dialog(tab = "Cooling",group="Geometry",enable = enable_setting_cooling_geometry));
  parameter .Modelica.Units.SI.Area A_cooling = pi*D_cooling*D_cooling/4 "Cross-section area of single cooling channel)" annotation (Dialog(tab = "Cooling",group="Geometry",enable = enable_setting_cooling_geometry));
  parameter .Modelica.Units.SI.Length C_cooling = pi*D_cooling "Circumference of single cooling channel" annotation (Dialog(tab = "Cooling",group="Geometry",enable = enable_setting_cooling_geometry));
  parameter .Modelica.Units.SI.Volume V_cooling = n_channels_cooling*A_cooling*L_cooling "Total volume of all cooling channels" annotation (Dialog(tab = "Cooling",group="Geometry",enable = enable_setting_cooling_geometry));
  parameter .Modelica.Units.SI.Area A_heat_cooling = C_cooling*L_cooling "Heat transfer area of single cooling channel" annotation (Dialog(tab = "Cooling",group="Geometry",enable = enable_setting_cooling_geometry));

  final parameter .Modelica.Units.SI.Length Dhyd_cooling = .Modelon.ThermoFluid.Functions.CharacteristicNumbers.HydraulicDiameter(A_cooling, C_cooling) "Hydraulic diameter of single cooling channel)" annotation (Dialog(group="Cooling channel geometry"));

  parameter Boolean includeStaticHead_cooling = false "If true, static head considered" annotation(Evaluate,Dialog(tab = "Cooling",group="Static head",enable = enable_setting_cooling_geometry));
  parameter .Modelica.Units.SI.Length height_cooling = 0 "Height(outlet) - height(inlet)" annotation (Dialog(enable=includeStaticHead_cooling, tab = "Cooling",group="Static head"));
  constant .Modelica.Units.SI.Acceleration g_cooling = .Modelica.Constants.g_n "Gravitational acceleration"annotation (Dialog(enable=false, tab = "Cooling",group="Static head"));

equation

  assert(L_cooling >= abs(height_cooling), "The parameter L must be >= abs(height).");

  annotation (Documentation(revisions="<html>
Copyright &copy; 2004-2024, MODELON AB <br /> The use of this software component is regulated by the licensing conditions for Modelon Libraries. <br />This copyright notice must, unaltered, accompany all components that are derived from, copied from, <br />or by other means have their origin from any Modelon Library.
</html>", info="<html>
<h4>LumpedCoolingGeometry</h4>
<p>Interface with geometry parameters for cooling channel.</p>
</html>"));

end LumpedCoolingGeometry;