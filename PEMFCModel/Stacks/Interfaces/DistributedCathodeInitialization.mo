within PEMFCModel.Stacks.Interfaces;
partial model DistributedCathodeInitialization "Initial parameters for distributed cathode channel"
  
  parameter Integer n(min=1) = 1 "Number of control volumes" annotation(Evaluate=true,Dialog(enable = true));
  parameter Boolean enable_setting_cathode_init = true "If true, selectability of cathode initialization parameter enabled" annotation (Dialog(tab="Initialization", group="Cathode channel"));
  parameter .Modelon.ThermoFluid.Choices.InitOptions initOpt_cathode = .Modelon.ThermoFluid.Choices.InitOptions.initialValues "Initialization option (steady state or fix value)" annotation(Dialog(tab="Initialization", group="Cathode channel",enable = enable_setting_cathode_init));
  parameter .Modelica.Units.SI.Pressure p_start_in_cathode = 1e5 "Initial pressure at cathode channel inlet" annotation (Dialog(tab="Initialization", group="Cathode channel",enable = enable_setting_cathode_init));
  parameter .Modelica.Units.SI.Pressure p_start_out_cathode = 1e5 "Initial pressure at cathode channel outlet" annotation (Dialog(tab="Initialization", group="Cathode channel",enable = enable_setting_cathode_init));
  parameter .Modelica.Units.SI.Pressure[n + 1] p_start_cathode=linspace(
      p_start_in_cathode,
      p_start_out_cathode,
      n + 1) "Initial pressure in control volumes and outlet port at cathode" annotation (Dialog(tab="Initialization", group="Cathode channel",enable = enable_setting_cathode_init));

  parameter Boolean initFromEnthalpy_cathode = false "If true, initialization from enthalpy, otherwise from temperature" annotation(Dialog(tab="Initialization", group="Cathode channel",enable = enable_setting_cathode_init));
  parameter .Modelica.Units.SI.SpecificEnthalpy h_start_in_cathode = 300e3 "Initial specific enthalpy at cathode channel inlet" annotation (Dialog(
      enable=initFromEnthalpy_cathode,
      tab="Initialization",
      group="Cathode channel"));
  parameter .Modelica.Units.SI.SpecificEnthalpy h_start_out_cathode = 300e3 "Initial specific enthalpy at cathode channel outlet" annotation (Dialog(
      enable=initFromEnthalpy_cathode,
      tab="Initialization",
      group="Cathode channel"));
  parameter .Modelica.Units.SI.SpecificEnthalpy[n + 1] h_start_cathode = linspace(
      h_start_in_cathode,
      h_start_out_cathode,
      n + 1) "Initial specific enthalpy in control volumes and inlet port at cathode" annotation (Dialog(
      enable=initFromEnthalpy_cathode,
      tab="Initialization",
      group="Cathode channel"));

  parameter .Modelica.Units.SI.Temperature T_start_in_cathode = 298.15 "Initial temperature at cathode channel inlet" annotation (
      Dialog(
      enable=not initFromEnthalpy_cathode,
      tab="Initialization",
      group="Cathode channel"));
  parameter .Modelica.Units.SI.Temperature T_start_out_cathode = 298.15 "Initial temperature at cathode channel outlet" annotation (
      Dialog(
      enable=not initFromEnthalpy_cathode,
      tab="Initialization",
      group="Cathode channel"));
  parameter .Modelica.Units.SI.Temperature[n + 1] T_start_cathode = linspace(
      T_start_in_cathode,
      T_start_out_cathode,
      n + 1) "Initial temperature in control volumes and inlet port at cathode" annotation (Dialog(
      enable=not initFromEnthalpy_cathode,
      tab="Initialization",
      group="Cathode channel"));

  parameter .Modelica.Units.SI.MassFraction X_start_cathode[:] "Initial mass-based composition in all control volumes at cathode" annotation (Dialog(tab="Initialization", group="Cathode channel"));

  parameter .Modelica.Units.SI.MassFlowRate m_flow_start_cathode = 0.1 "Initial mass flow rate (guessed value) at cathode" annotation (Dialog(tab="Initialization", group="Cathode channel",enable = enable_setting_cathode_init));

  annotation (Documentation(revisions="<html>
Copyright &copy; 2004-2024, MODELON AB <br /> The use of this software component is regulated by the licensing conditions for Modelon Libraries. <br />This copyright notice must, unaltered, accompany all components that are derived from, copied from, <br />or by other means have their origin from any Modelon Library.
</html>", info="<html>
<h4>DistributedCathodeInitialization</h4>
<p>Interface with cathode initialialization parameters.</p>
</html>"));

end DistributedCathodeInitialization;
