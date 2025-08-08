within PEMFCModel.Membrane.Losses.ConcentrationLoss;

model ZeroLoss "Zero concentration loss"
  
  extends PEMFCModel.Membrane.Losses.Templates.BaseConcentrationLoss;

equation
    
  for i in 1:N loop
    
    E_loss_cell[i] = 0;
    
  end for;

  annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
                             Text(
          extent={{-80,20},{80,-20}},
          textColor={28,108,200},
          textString="ZeroLoss")}),                              Diagram(coordinateSystem(preserveAspectRatio=false)),
    Documentation(revisions="<html>
Copyright &copy; 2004-2024, MODELON AB <br /> The use of this software component is regulated by the licensing conditions for Modelon Libraries. <br />This copyright notice must, unaltered, accompany all components that are derived from, copied from, <br />or by other means have their origin from any Modelon Library.
</html>"));

end ZeroLoss;
