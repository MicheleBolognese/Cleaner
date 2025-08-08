within PEMFCModel.Membrane.MassTransport.WaterPermeation;
model ZeroFlow "Zero flow"
  
  extends PEMFCModel.Membrane.MassTransport.Templates.PartialWaterPermeation(m_flow = zeros(N));

equation
  
  k = zeros(N);

  annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)),
    Documentation(info="<html>
<p>Sets the Darcy water permeation coefficient and mass flow rate to zero.</p>
</html>", revisions="<html>
Copyright &copy; 2004-2024, MODELON AB <br /> The use of this software component is regulated by the licensing conditions for Modelon Libraries. <br />This copyright notice must, unaltered, accompany all components that are derived from, copied from, <br />or by other means have their origin from any Modelon Library.
</html>"));
    
end ZeroFlow;
