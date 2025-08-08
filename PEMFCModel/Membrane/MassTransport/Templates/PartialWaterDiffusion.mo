within PEMFCModel.Membrane.MassTransport.Templates;
partial model PartialWaterDiffusion

  extends PEMFCModel.Membrane.MassTransport.Interfaces.PartialWaterFlow(m_flow = {-D_diff[i]*(c_cath_H2O[i] - c_an_H2O[i])/z*A_cell/N*MM_H2O*n_cell for i in 1:N});

  input Real lambda_an[N] "Membrane water content at the anode-side in each node" annotation(Dialog(enable = enableInternal, tab = "Internal", group = "Input variables"));
  input Real lambda_cath[N] "Membrane water content at the cathode-side in each node" annotation(Dialog(enable = enableInternal, tab = "Internal", group = "Input variables"));
  
  .Modelica.Units.SI.DiffusionCoefficient D_diff[N] "Water diffusion coefficient through membrane in each node";
  .Modelica.Units.SI.Concentration c_an_H2O[N] "H2O concentration at the anode side in each node";
  .Modelica.Units.SI.Concentration c_cath_H2O[N] "H2O concentration at the cathode side in each node";

  annotation (
    defaultComponentName="membraneWaterProps",
    Icon(coordinateSystem(preserveAspectRatio=false), graphics={Rectangle(
          extent={{-100,100},{100,-100}},
          lineColor={0,0,0},
          fillPattern=FillPattern.Solid,
          fillColor={255,255,255}), Text(
          extent={{-40,20},{40,-20}},
          textColor={0,0,0},
          textString="WD")}),
    Diagram(coordinateSystem(preserveAspectRatio=false)),
    Documentation(revisions="<html>
Copyright &copy; 2004-2024, MODELON AB <br /> The use of this software component is regulated by the licensing conditions for Modelon Libraries. <br />This copyright notice must, unaltered, accompany all components that are derived from, copied from, <br />or by other means have their origin from any Modelon Library.
</html>", info="<html>
<p>Template class for the definition of diffusion of water across a membrane. The default equation is based on the concentration difference across the membrane (i.e., Fick&apos;s law).</p>
</html>"));

end PartialWaterDiffusion;
