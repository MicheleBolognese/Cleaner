within PEMFCModel.Membrane.Losses.Interfaces;

partial model LossInterface "Two-pin interface of any membrane voltage loss component"

  parameter Boolean enableInternal = true annotation(Dialog(tab = "Internal", group = "Input parameters"));

  parameter Integer n_cell(min = 1) = 50 "Number of series-connected cells" annotation(Dialog(enable = enableInternal, tab = "Internal", group = "Input parameters"));
  parameter Integer N(min = 1) = 1 "Number of along-the-channel discretization nodes" annotation(Dialog(enable = enableInternal, tab = "Internal", group = "Input parameters"));

  .Modelica.Electrical.Analog.Interfaces.PositivePin pin_p[N] annotation (Placement(transformation(extent={{-110,-10},{-90,10}})));
  .Modelica.Electrical.Analog.Interfaces.NegativePin pin_n[N] annotation (Placement(transformation(extent={{90,-10},{110,10}})));
  
  annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={
          Rectangle(extent={{-100,40},{100,-40}}, lineColor={28,108,200})}),
      Diagram(coordinateSystem(preserveAspectRatio=false)),
    Documentation(info="<html>
Two-pin interface for membrane loss models.
</html>", revisions="<html>
Copyright &copy; 2004-2024, MODELON AB <br /> The use of this software component is regulated by the licensing conditions for Modelon Libraries. <br />This copyright notice must, unaltered, accompany all components that are derived from, copied from, <br />or by other means have their origin from any Modelon Library.
</html>"));

end LossInterface;
