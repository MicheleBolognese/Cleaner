within PEMFCModel.Membrane.Losses.Templates;
partial model BaseConcentrationLoss "Template for concentration loss"
  
  extends PEMFCModel.Membrane.Losses.Templates.BaseLoss;
  
  input .Modelica.Units.SI.CurrentDensity j_ionic[N] "Current density in each node" annotation(Dialog(enable = enableInternal, tab = "Internal", group = "Input parameters"));

  annotation (Documentation(revisions="<html>
Copyright &copy; 2004-2024, MODELON AB <br /> The use of this software component is regulated by the licensing conditions for Modelon Libraries. <br />This copyright notice must, unaltered, accompany all components that are derived from, copied from, <br />or by other means have their origin from any Modelon Library.
</html>"));

end BaseConcentrationLoss;