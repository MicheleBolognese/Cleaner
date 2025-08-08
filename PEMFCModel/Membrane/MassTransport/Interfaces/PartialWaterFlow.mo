within PEMFCModel.Membrane.MassTransport.Interfaces;
partial model PartialWaterFlow

  extends PEMFCModel.Membrane.MassTransport.Interfaces.PartialWaterModel;
  
  input Real lambda[N] "Membrane water content" annotation(Dialog(enable = enableInternal, tab = "Internal", group = "Input variables"));
  
  .Modelica.Units.SI.MassFlowRate m_flow[N] = zeros(N) "Mass flow rate through the membrane in each node (> 0 from anode to cathode)";

  annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)),
    Documentation(info="<html>
Interface containing additional parameters/variables for water transport model.
</html>", revisions="<html>
Copyright &copy; 2004-2024, MODELON AB <br /> The use of this software component is regulated by the licensing conditions for Modelon Libraries. <br />This copyright notice must, unaltered, accompany all components that are derived from, copied from, <br />or by other means have their origin from any Modelon Library.
</html>"));

end PartialWaterFlow;
