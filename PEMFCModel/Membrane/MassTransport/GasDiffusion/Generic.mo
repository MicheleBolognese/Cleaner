within PEMFCModel.Membrane.MassTransport.GasDiffusion;

model Generic "User defined"
  
  extends PEMFCModel.Membrane.MassTransport.Templates.PartialGasDiffusion;

  input Real K_perm0[N, nDiff](unit = "mol/(s.m.Pa)") = fill(1e-14, N, nDiff) "Permeability coefficient of ech diffusive species in each node (diffusion*solubility)" annotation (Dialog(group = "Input variables"));

equation

  K_perm = CF_N.*K_perm0;

  annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)),
    Documentation(info="<html>
<p>The simplest definition for the permeability coefficient. Allows for a custom formulation to be specified either directly or referenced from some other component for each gas species added. Input may be constant or time-dependent.</p>
</html>", revisions="<html>
Copyright &copy; 2004-2024, MODELON AB <br /> The use of this software component is regulated by the licensing conditions for Modelon Libraries. <br />This copyright notice must, unaltered, accompany all components that are derived from, copied from, <br />or by other means have their origin from any Modelon Library.
</html>"));

end Generic;
