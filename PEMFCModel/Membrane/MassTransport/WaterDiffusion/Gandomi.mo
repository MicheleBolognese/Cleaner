within PEMFCModel.Membrane.MassTransport.WaterDiffusion;
model Gandomi "Correlation by Gandomi for water diffusion"

  extends PEMFCModel.Membrane.MassTransport.Templates.PartialWaterDiffusion(enableInternal = false);

protected
  
  Real lambda_int[N];

equation

  for i in 1:N loop
    
    c_an_H2O[i] = p_an_H2O[i]/(.Modelica.Constants.R*T_an[i]);
    c_cath_H2O[i] = p_cath_H2O[i]/(.Modelica.Constants.R*T_cath[i]);
    
    lambda_int[i] = .Modelon.Math.Smoothing.smoothBounds(lambda[i], {0,17}, {0.01,16.9});

    D_diff[i] = CF_N[i]*.Modelon.Math.Smoothing.spliceFunction(
      4.17e-8*lambda_int[i]*(1 + 161*exp(-lambda_int[i]))*exp(-2436/T_cell[i]),
      3.10e-7*lambda_int[i]*(-1 + exp(0.28*lambda_int[i]))*exp(-2436/T_cell[i]),
      lambda_int[i] - 3, 0.5);
  
  end for;

  annotation (
    defaultComponentName="membraneWaterProps",
    Icon(coordinateSystem(preserveAspectRatio=false)),
    Diagram(coordinateSystem(preserveAspectRatio=false)),
    Documentation(info="<html>
<p>Diffusion coefficient approach which is a function of membrane water content as defined by Gandomi et. al.</p>
<h4>References</h4>
<p>Gandomi, Yasser Ashraf, M. D. Edmundson, F. C. Busby, and Matthew M. Mench. &ldquo;Water Management in Polymer Electrolyte Fuel Cells through Asymmetric Thermal and Mass Transport Engineering of the Micro-Porous Layers.&rdquo; <i>Journal of The Electrochemical Society</i> 163, no. 8 (June 16, 2016): F933. <a href=\"https://doi.org/10.1149/2.1331608jes\">https://doi.org/10.1149/2.1331608jes</a>. </p>
</html>", revisions="<html>
Copyright &copy; 2004-2024, MODELON AB <br /> The use of this software component is regulated by the licensing conditions for Modelon Libraries. <br />This copyright notice must, unaltered, accompany all components that are derived from, copied from, <br />or by other means have their origin from any Modelon Library.
</html>"));

end Gandomi;
