within PEMFCModel.Stacks.Interfaces;
partial model LumpedAnodeGeometry "Anode channel geometry parameters"

  constant Real pi = .Modelica.Constants.pi;
  parameter Boolean enable_setting_anode_geometry = true "If true, selectability of anode geometry parameters enabled" annotation(Evaluate=true,Dialog(tab = "Substack",group="Anode side geometry"));

  parameter Real n_channels_anode(min=1.0) = 1.0 "Number of parallel channels at anode" annotation(Dialog(tab = "Substack", group="Anode side geometry",enable = enable_setting_anode_geometry));
  parameter .Modelica.Units.SI.Length L_anode = 0.1 "Length of anode channel" annotation (Dialog(tab = "Substack",group="Anode side geometry",enable = enable_setting_anode_geometry));
  parameter .Modelica.Units.SI.Length D_anode = 0.01 "Diameter of single anode channel" annotation (Dialog(tab = "Substack",group="Anode side geometry",enable = enable_setting_anode_geometry));
  parameter .Modelica.Units.SI.Length C_anode = pi*D_anode "Circumference of single anode channel" annotation(Dialog(tab = "Substack", group="Anode side geometry",enable = enable_setting_anode_geometry));
  parameter .Modelica.Units.SI.Area A_anode=pi*D_anode*D_anode/4 "Cross section area of single anode channel" annotation(Dialog(tab = "Substack", group="Anode side geometry",enable = enable_setting_anode_geometry));
  parameter .Modelica.Units.SI.Volume V_anode = n_channels_anode*L_anode*A_anode "Total volume of anode channels" annotation(Dialog(tab = "Substack", group="Anode side geometry",enable = enable_setting_anode_geometry));
  parameter .Modelica.Units.SI.Area A_heat_anode = C_anode*L_anode "Heat transfer area of single anode channel" annotation(Dialog(tab = "Substack", group="Anode side geometry",enable = enable_setting_anode_geometry));
  
  final parameter .Modelica.Units.SI.Length Dhyd_anode = .Modelon.ThermoFluid.Functions.CharacteristicNumbers.HydraulicDiameter(A_anode, C_anode) "Hydraulic diameter of single channel" annotation (Dialog(group="Anode side geometry",enable = true));

  annotation (Documentation(revisions="<html>
Copyright &copy; 2004-2024, MODELON AB <br /> The use of this software component is regulated by the licensing conditions for Modelon Libraries. <br />This copyright notice must, unaltered, accompany all components that are derived from, copied from, <br />or by other means have their origin from any Modelon Library.
</html>", info="<html>
<h4>LumpedAnodeGeometry</h4>
<p>Interface with anode channel geometry parameters.</p>
</html>"));

end LumpedAnodeGeometry;