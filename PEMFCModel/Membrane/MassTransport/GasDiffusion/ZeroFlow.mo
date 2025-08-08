within PEMFCModel.Membrane.MassTransport.GasDiffusion;
model ZeroFlow "Zero flow"
  
  extends PEMFCModel.Membrane.MassTransport.Templates.PartialGasDiffusion(mX_flow_diff = zeros(N, nDiff));

equation

  K_perm = zeros(N, nDiff);

  annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)),
    Documentation(info="<html>
<p>Sets the permeability coefficient to zero which yields in no transport for all gas species.</p>
</html>", revisions="<html>
Copyright &copy; 2004-2024, MODELON AB <br /> The use of this software component is regulated by the licensing conditions for Modelon Libraries. <br />This copyright notice must, unaltered, accompany all components that are derived from, copied from, <br />or by other means have their origin from any Modelon Library.
</html>"));

end ZeroFlow;
