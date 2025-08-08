within PEMFCModel.Pipes.Interfaces;
partial model DistributedChannelInitialization "Initialization parameters for distributed pipe"

  parameter Integer ni(min=1) = 1 "Number of control volumes" annotation(Evaluate=true,Dialog(group="Discretization"));

/* Initialization */
  parameter .Modelon.ThermoFluid.Choices.InitOptions initOpt = .Modelon.ThermoFluid.Choices.InitOptions.initialValues "Initialization option (steady state or fix value)" annotation(Dialog(tab="Initialization", group="Initialization option"));
  parameter .Modelica.Units.SI.Pressure p_start_in = 1e5 "Initial inlet pressure" annotation (Dialog(tab="Initialization", group="Pressure"));
  parameter .Modelica.Units.SI.Pressure p_start_out = 1e5 "Initial outlet pressure" annotation (Dialog(tab="Initialization", group="Pressure"));
  parameter .Modelica.Units.SI.Pressure p_start[ni + 1] = linspace(
      p_start_in,
      p_start_out,
      ni + 1) "Initial pressure in control volumes and outlet port" annotation (Dialog(tab="Initialization", group="Pressure"));

  parameter Boolean initFromEnthalpy = false "If true, initialization from enthalpy, otherwise from temperature" annotation(Dialog(tab="Initialization", group="Enthalpy"), Evaluate=true);
  parameter .Modelica.Units.SI.SpecificEnthalpy h_start_in = 300e3 "Initial inlet specific enthalpy" annotation (Dialog(
      enable=initFromEnthalpy,
      tab="Initialization",
      group="Enthalpy"),
    __Modelon(DialogExtensions(hide= not initFromEnthalpy)));
  parameter .Modelica.Units.SI.SpecificEnthalpy h_start_out = 300e3 "Initial outlet specific enthalpy" annotation (Dialog(
      enable=initFromEnthalpy,
      tab="Initialization",
      group="Enthalpy"),
    __Modelon(DialogExtensions(hide= not initFromEnthalpy)));
  parameter .Modelica.Units.SI.SpecificEnthalpy[ni + 1] h_start = linspace(
      h_start_in,
      h_start_out,
      ni + 1) "Initial specific enthalpy in control volumes and inlet port" annotation (Dialog(
      enable=initFromEnthalpy,
      tab="Initialization",
      group="Enthalpy"));

  parameter .Modelon.Media.Units.Temperature T_start_in = 298.15 "Initial inlet temperature" annotation(Dialog(enable=not initFromEnthalpy, tab="Initialization", group="Temperature"),
    __Modelon(DialogExtensions(hide= initFromEnthalpy)));
  parameter .Modelon.Media.Units.Temperature T_start_out = 298.15 "Initial outlet temperature" annotation(Dialog(enable=not initFromEnthalpy, tab="Initialization", group="Temperature"),
    __Modelon(DialogExtensions(hide= initFromEnthalpy)));
  parameter .Modelica.Units.SI.Temperature T_start[ni + 1] = linspace(
      T_start_in,
      T_start_out,
      ni + 1) "Initial temperature in control volumes and inlet port" annotation (Dialog(
      enable=not initFromEnthalpy,
      tab="Initialization",
      group="Temperature"));

  parameter .Modelica.Units.SI.MassFraction X_start[:] "Initial mass-based composition in all control volumes" annotation (Dialog(tab="Initialization", group="Mass fractions"));

  parameter Real C_start[:] "Initial trace component concentrations in all control volumes" annotation(Dialog(tab="Initialization", group="Trace component fractions"));

  parameter .Modelica.Units.SI.MassFlowRate m_flow_start = 0.1 "Initial mass flow rate (guess value)"annotation (Dialog(tab="Initialization", group="Mass flow rate"));

end DistributedChannelInitialization;