within PEMFCModel.Membrane.MassTransport.Interfaces;
partial model PartialTransport

  parameter Boolean enableInternal = true " = true to enable selectability of 'Internal' parameters" annotation(Dialog(tab = "Internal", group = "Input parameters"));
  
  replaceable package Medium_an = .FuelCell.Media.PreDefined.IdealGases.NASAReformate constrainedby .FuelCell.Media.Templates.ReactionGas "Anodic medium" annotation(choicesAllMatching, Dialog(tab = "Internal", enable = enableInternal, group = "Input parameters"));
  replaceable package Medium_cath = .FuelCell.Media.PreDefined.IdealGases.NASAMoistAir constrainedby .FuelCell.Media.Templates.ReactionGas "Cathodic medium" annotation(choicesAllMatching, Dialog(tab = "Internal", enable = enableInternal, group = "Input parameters"));

  parameter Integer N(min = 1) = 1 "Number of along-the-channel discretization nodes" annotation(Dialog(tab = "Internal", enable = enableInternal, group = "Input parameters"));
  parameter Integer n_cell(min = 1) = 1 "Number of series-connected cells" annotation(Dialog(tab = "Internal", enable = enableInternal, group = "Input parameters"));
  parameter .Modelica.Units.SI.Area A_cell = 180e-4 "Cell active area" annotation(Dialog(tab = "Internal", enable = enableInternal, group = "Input parameters"));
  parameter .Modelica.Units.SI.Length z = 20e-6 "Membrane thickness" annotation(Dialog(tab = "Internal", enable = enableInternal, group = "Input parameters"));
  parameter .Modelica.Units.SI.Density rho_dry_m = 1980 "Dry membrane density" annotation(Dialog(tab = "Internal", enable = enableInternal, group = "Input parameters"));
  parameter .Modelica.Units.SI.MolarMass EW_m = 1.1 "Membrane equivalent weight" annotation(Dialog(tab = "Internal", enable = enableInternal, group = "Input parameters"));

  input .Modelica.Units.SI.Pressure p_an_partial[N,Medium_an.nS] "Partial pressure of each chemical species at the anode in each node" annotation(Dialog(tab = "Internal", enable = enableInternal, group = "Input variables"));
  input .Modelica.Units.SI.Pressure p_cath_partial[N,Medium_cath.nS] "Partial pressure of each chemical species at the cathode in each node" annotation(Dialog(tab = "Internal", enable = enableInternal, group = "Input variables"));
  input .Modelica.Units.SI.Temperature T_an[N] "Temperature at the anode side in each node" annotation(Dialog(tab = "Internal", enable = enableInternal, group = "Input variables"));
  input .Modelica.Units.SI.Temperature T_cath[N] "Temperature at the cathode side in each node" annotation(Dialog(tab = "Internal", enable = enableInternal, group = "Input variables"));
  input .Modelica.Units.SI.Temperature T_cell[N] "Temperature of each cell node" annotation(Dialog(tab = "Internal", enable = enableInternal, group = "Input variables"));

  annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)),
    Documentation(info="<html>
Interface containing parameters and variables declarations for mass transport models.
</html>", revisions="<html>
Copyright &copy; 2004-2024, MODELON AB <br /> The use of this software component is regulated by the licensing conditions for Modelon Libraries. <br />This copyright notice must, unaltered, accompany all components that are derived from, copied from, <br />or by other means have their origin from any Modelon Library.
</html>"));

end PartialTransport;
