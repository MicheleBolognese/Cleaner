within PEMFCModel.Stacks.Interfaces;
partial model CoolingInitialization "Initialization parameters for distributed cooling channel"

  parameter Integer ni_cooling(min = 1) = 1 "Number of control volumes" annotation(Evaluate=true,Dialog(enable = true));
  
  parameter Boolean enable_setting_cooling_init = true "If true, selectability of cooling channel initialization parameters enabled" annotation(Dialog(tab = "Initialization",group = "Cooling channel"));
/* Initialization */
  parameter .Modelon.ThermoFluid.Choices.InitOptions initOpt_cooling = .Modelon.ThermoFluid.Choices.InitOptions.initialValues "Initialization option (steady state or fix value)" annotation(Dialog(tab="Initialization", group="Cooling channel",enable = enable_setting_cooling_init));
  parameter .Modelica.Units.SI.Pressure p_start_in_cooling = 1e5 "Initial inlet pressure" annotation (Dialog(tab="Initialization", group="Cooling channel",enable = enable_setting_cooling_init));
  parameter .Modelica.Units.SI.Pressure p_start_out_cooling = 1e5 "Initial outlet pressure" annotation (Dialog(tab="Initialization", group="Cooling channel",enable = enable_setting_cooling_init));
  parameter .Modelica.Units.SI.Pressure p_start_cooling[ni_cooling + 1] = linspace(
    p_start_in_cooling,
    p_start_out_cooling,
    ni_cooling + 1) "Initial pressure in control volumes and outlet port" annotation (Dialog(tab="Initialization", group="Cooling channel",enable = enable_setting_cooling_init));

  parameter Boolean initFromEnthalpy_cooling = false "If true, initialization from enthalpy, otherwise from temperature" annotation(Dialog(tab="Initialization", group="Cooling channel",enable = enable_setting_cooling_init));
  parameter .Modelica.Units.SI.SpecificEnthalpy h_start_in_cooling = 300e3 "Initial inlet specific enthalpy" annotation (Dialog(     enable=initFromEnthalpy_cooling,tab="Initialization", group="Cooling channel"));
  parameter .Modelica.Units.SI.SpecificEnthalpy h_start_out_cooling = 300e3 "Initial outlet specific enthalpy" annotation (Dialog(
      enable=initFromEnthalpy_cooling,
      tab="Initialization", group="Cooling channel"));
  parameter .Modelica.Units.SI.SpecificEnthalpy h_start_cooling[ni_cooling + 1] = linspace(
    h_start_in_cooling,
    h_start_out_cooling,
    ni_cooling + 1) "Initial specific enthalpy in  control volumes and inlet port" annotation (Dialog(
      enable=initFromEnthalpy_cooling,
      tab="Initialization", group="Cooling channel"));

  parameter .Modelica.Units.SI.Temperature T_start_in_cooling = 298.15 "Initial inlet temperature" annotation (Dialog(
      enable=not initFromEnthalpy_cooling,
      tab="Initialization", group="Cooling channel"));
  parameter .Modelica.Units.SI.Temperature T_start_out_cooling=298.15 "Initial outlet temperature" annotation (Dialog(
      enable=not initFromEnthalpy_cooling,
      tab="Initialization", group="Cooling channel"));
  parameter .Modelica.Units.SI.Temperature T_start_cooling[ni_cooling + 1] = linspace(
    T_start_in_cooling,
    T_start_out_cooling,
    ni_cooling + 1) "Initial temperature in  control volumes and inlet port" annotation (Dialog(
      enable=not initFromEnthalpy_cooling,
      tab="Initialization", group="Cooling channel"));

  parameter .Modelica.Units.SI.MassFraction[:] X_start_cooling "Initial mass-based composition in all control volumes " annotation (Dialog( tab="Initialization", group="Cooling channel"));

  parameter .Modelica.Units.SI.MassFlowRate m_flow_start_cooling = 0.1 "Initial mass flow rate (guess value)" annotation (Dialog(tab="Initialization", group="Cooling channel",enable = enable_setting_cooling_init));

  annotation (Documentation(revisions="<html>
Copyright &copy; 2004-2024, MODELON AB <br /> The use of this software component is regulated by the licensing conditions for Modelon Libraries. <br />This copyright notice must, unaltered, accompany all components that are derived from, copied from, <br />or by other means have their origin from any Modelon Library.
</html>", info="<html>
<h4>CoolingInitParameters</h4>
<p>Interface with initialization parameters for cooling pipe.</p>
</html>"));

end CoolingInitialization;