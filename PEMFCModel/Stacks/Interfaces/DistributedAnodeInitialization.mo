within PEMFCModel.Stacks.Interfaces;
partial model DistributedAnodeInitialization "Initial parameters for distributed anode channel"
  
  parameter Integer n(min = 1) = 1 "Number of control volumes" annotation(Evaluate=true,Dialog(enable = true));
  parameter Boolean enable_setting_anode_init = true "If true, selectability of anode initialization parameters enabled" annotation(Dialog(tab="Initialization", group="Anode channel"));
  parameter .Modelon.ThermoFluid.Choices.InitOptions initOpt_anode = .Modelon.ThermoFluid.Choices.InitOptions.initialValues "Initialization option (steady state or fix value)" annotation(Dialog(tab="Initialization", group="Anode channel",enable = enable_setting_anode_init));
  parameter .Modelica.Units.SI.Pressure p_start_in_anode = 1e5 "Initial pressure at anode channel inlet" annotation (Dialog(tab="Initialization", group="Anode channel",enable = enable_setting_anode_init));
  parameter .Modelica.Units.SI.Pressure p_start_out_anode = 1e5 "Initial pressure at anode channel outlet " annotation (Dialog(tab="Initialization", group="Anode channel",enable = enable_setting_anode_init));
  parameter .Modelica.Units.SI.Pressure[n + 1] p_start_anode = linspace(
      p_start_in_anode,
      p_start_out_anode,
      n + 1) "Initial pressure in control volumes and outlet port at anode" annotation (Dialog(tab="Initialization", group="Anode channel",enable = enable_setting_anode_init));

  parameter Boolean initFromEnthalpy_anode = false "If true, initialization from enthalpy, otherwise from temperature" annotation(Dialog(tab="Initialization", group="Anode channel",enable = enable_setting_anode_init));
  parameter .Modelica.Units.SI.SpecificEnthalpy h_start_in_anode = 300e3 "Initial specific enthalpy at anode channel inlet" annotation (Dialog(
      enable=initFromEnthalpy_anode,
      tab="Initialization",
      group="Anode channel"));
  parameter .Modelica.Units.SI.SpecificEnthalpy h_start_out_anode = 300e3 "Initial specific enthalpy at anode channel outlet" annotation (Dialog(
      enable=initFromEnthalpy_anode,
      tab="Initialization",
      group="Anode channel"));
  parameter .Modelica.Units.SI.SpecificEnthalpy[n + 1] h_start_anode = linspace(
      h_start_in_anode,
      h_start_out_anode,
      n + 1) "Initial specific enthalpy in control volumes and inlet port at anode" annotation (Dialog(
      enable=initFromEnthalpy_anode,
      tab="Initialization",
      group="Anode channel"));

  parameter .Modelica.Units.SI.Temperature T_start_in_anode = 298.15 "Initial temperature at anode channel inlet" annotation (Dialog(
      enable=not initFromEnthalpy_anode,
      tab="Initialization",
      group="Anode channel"));
  parameter .Modelica.Units.SI.Temperature T_start_out_anode = 298.15 "Initial temperature at anode channel outlet" annotation (
      Dialog(
      enable=not initFromEnthalpy_anode,
      tab="Initialization",
      group="Anode channel"));
  parameter .Modelica.Units.SI.Temperature[n + 1] T_start_anode = linspace(
      T_start_in_anode,
      T_start_out_anode,
      n + 1) "Initial temperature in control volumes and inlet port at anode" annotation (Dialog(
      enable=not initFromEnthalpy_anode,
      tab="Initialization",
      group="Anode channel"));

  parameter .Modelica.Units.SI.MassFraction X_start_anode[:] "Initial mass-based composition in all control volumes at anode" annotation (Dialog(tab="Initialization", group="Anode channel"));

  parameter .Modelica.Units.SI.MassFlowRate m_flow_start_anode = 0.1 "Initial mass flow rate (guessed value) at anode" annotation (Dialog(tab="Initialization", group="Anode channel",enable = enable_setting_anode_init));

  annotation (Documentation(revisions="<html>
Copyright &copy; 2004-2024, MODELON AB <br /> The use of this software component is regulated by the licensing conditions for Modelon Libraries. <br />This copyright notice must, unaltered, accompany all components that are derived from, copied from, <br />or by other means have their origin from any Modelon Library.
</html>", info="<html>
<h4>DistributedAnodeInitialization</h4>
<p>Interface with anode initialization parameters.</p>
</html>"));

end DistributedAnodeInitialization;