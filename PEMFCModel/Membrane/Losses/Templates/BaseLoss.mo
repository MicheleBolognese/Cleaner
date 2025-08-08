within PEMFCModel.Membrane.Losses.Templates;
partial model BaseLoss "Base class for fuel cell voltage loss"
  
  extends PEMFCModel.Membrane.Losses.Interfaces.LossInterface;

  .Modelica.Units.SI.Voltage E_loss_stack[N](each min = 0, start = 0) "Stack voltage loss in each node";
  .Modelica.Units.SI.Voltage E_loss_cell[N](each min = 0, start = 0) "Single cell voltage loss in each node";
  .Modelica.Units.SI.Power powerLoss_stack "Stack power loss";

equation

  pin_p.v - pin_n.v = E_loss_stack;
  fill(0.0, N) = pin_p.i + pin_n.i;
  E_loss_stack = E_loss_cell*n_cell;
  powerLoss_stack = sum(pin_p.i.*E_loss_stack)
  
  annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
        coordinateSystem(preserveAspectRatio=false)),
    Documentation(info="<html>
<p>Minimum template of membrane loss model. This can be used for adding user-defined loss model, by simply extending this template
and adding formulation for <code>E_loss</code>.</p> 
<p>It will then be multiplied by <code>n_cell</code> parameter to give the total voltage loss between the pins.</p>
</html>", revisions="<html>
Copyright &copy; 2004-2024, MODELON AB <br /> The use of this software component is regulated by the licensing conditions for Modelon Libraries. <br />This copyright notice must, unaltered, accompany all components that are derived from, copied from, <br />or by other means have their origin from any Modelon Library.
</html>"));

end BaseLoss;