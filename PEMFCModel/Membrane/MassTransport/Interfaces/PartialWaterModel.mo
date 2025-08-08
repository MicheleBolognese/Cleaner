within PEMFCModel.Membrane.MassTransport.Interfaces;
partial model PartialWaterModel

  extends PEMFCModel.Membrane.MassTransport.Interfaces.PartialTransport;
  
  package Medium_H2O = .Modelon.Media.PreDefined.TwoPhase.WaterIF97 "Medium for pure H2O properties" annotation(Dialog(tab = "Internal", enable = false, group = "Predefined: do not override"));
   
  parameter Integer iH2O_an = 1 "Index of H2O in anodic medium" annotation(Dialog(tab = "Internal", enable = enableInternal, group = "Input parameters")); 
  parameter Integer iH2O_cath = 1 "Index of H2O in cathodic medium" annotation(Dialog(tab = "Internal", enable = enableInternal, group = "Input parameters"));
  
  Real CF_0 = 1.0 "Correction factor (1 by default)" annotation(Dialog(tab = "Advanced", group = "Input parameters"));
  Real CF_N[N] = fill(CF_0, N) "Set if CF_0 is non-uniform" annotation(Dialog(tab="Advanced", group = "Input parameters"));

protected
  
  constant .Modelica.Units.SI.MolarMass MM_H2O = Medium_H2O.MMX[1] "H2O molar mass";
  .Modelica.Units.SI.Pressure p_an_H2O[N] = {p_an_partial[i, iH2O_an] for i in 1:N} "H2O partial pressure in anodic medium in each node";
  .Modelica.Units.SI.Pressure p_cath_H2O[N] = {p_cath_partial[i, iH2O_cath] for i in 1:N} "H2O partial pressure in anodic medium in each node";

  annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)),
    Documentation(info="<html>
Interface containing additional parameters/variables for water model.
</html>", revisions="<html>
Copyright &copy; 2004-2024, MODELON AB <br /> The use of this software component is regulated by the licensing conditions for Modelon Libraries. <br />This copyright notice must, unaltered, accompany all components that are derived from, copied from, <br />or by other means have their origin from any Modelon Library.
</html>"));

end PartialWaterModel;
