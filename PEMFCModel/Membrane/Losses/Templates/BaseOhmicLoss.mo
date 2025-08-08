within PEMFCModel.Membrane.Losses.Templates;
partial model BaseOhmicLoss "Template for ohmic loss"
  
  extends PEMFCModel.Membrane.Losses.Templates.BaseLoss;
  
  parameter .Modelica.Units.SI.Area A_cell = 180e-4 "Cell active area" annotation(Dialog(enable = enableInternal, tab = "Internal", group = "Input parameters"));
  
  .FuelCell.Internal.Units.AreaSpecificResistance ASR[N] "Area specific resistance in each node";
  .Modelica.Units.SI.Current I[N] "Current in each node";

equation
 
  for j in 1:N loop

    E_loss_cell[j] = ASR[j]*I[j]/(A_cell/N);
    I[j] = pin_p[j].i;

  end for;

  annotation (Documentation(revisions="<html>
Copyright &copy; 2004-2024, MODELON AB <br /> The use of this software component is regulated by the licensing conditions for Modelon Libraries. <br />This copyright notice must, unaltered, accompany all components that are derived from, copied from, <br />or by other means have their origin from any Modelon Library.
</html>"));

end BaseOhmicLoss;