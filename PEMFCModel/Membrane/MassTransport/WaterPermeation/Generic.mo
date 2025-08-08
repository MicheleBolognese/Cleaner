within PEMFCModel.Membrane.MassTransport.WaterPermeation;
model Generic "User defined"
  
  extends PEMFCModel.Membrane.MassTransport.Templates.PartialWaterPermeation;

  input .Modelica.Units.SI.Area k_0 = 1.58e-18 "Water Darcy's permeation coefficient" annotation(Dialog(group = "Input variables"));
  input .Modelica.Units.SI.Area k_N[N] = fill(k_0, N) "Set if k_0 is non-uniform" annotation(Dialog(group = "Input variables"));

equation
  
  k = CF_N.*k_N;

  annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)),
    Documentation(info="<html>
<p>The simplest definition for the darcy water permeation coefficient. Allows for a custom formulation to be specified either directly or referenced from some other component. Input may be constant or time-dependent.</p>
</html>", revisions="<html>
Copyright &copy; 2004-2024, MODELON AB <br /> The use of this software component is regulated by the licensing conditions for Modelon Libraries. <br />This copyright notice must, unaltered, accompany all components that are derived from, copied from, <br />or by other means have their origin from any Modelon Library.
</html>"));
    
end Generic;
