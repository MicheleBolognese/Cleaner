within PEMFCModel.Membrane.MassTransport.WaterDiffusion;
model ZeroFlow "Zero flow"

  extends PEMFCModel.Membrane.MassTransport.Templates.PartialWaterDiffusion(m_flow = zeros(N), c_an_H2O = zeros(N), c_cath_H2O = zeros(N));

equation

  D_diff = zeros(N);

  annotation (defaultComponentName="membraneWaterProps",Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)),
    Documentation(revisions="<html>
Copyright &copy; 2004-2024, MODELON AB <br /> The use of this software component is regulated by the licensing conditions for Modelon Libraries. <br />This copyright notice must, unaltered, accompany all components that are derived from, copied from, <br />or by other means have their origin from any Modelon Library.
</html>", info="<html>
<p>Sets the water diffusion coefficient and mass flow rate to zero.</p>
</html>"));
    
end ZeroFlow;
