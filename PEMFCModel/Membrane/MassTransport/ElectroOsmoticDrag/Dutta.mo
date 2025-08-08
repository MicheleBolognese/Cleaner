within PEMFCModel.Membrane.MassTransport.ElectroOsmoticDrag;
model Dutta "Correlation by Dutta for electro-osmotic drag"
  
  extends PEMFCModel.Membrane.MassTransport.Templates.PartialElectroOsmoticDrag(enableInternal = false);

  final parameter Real coeffs[3] = {-3.4e-19,0.05,0.0029} "Coefficients for c[1] + c[2]*lambda + ..." annotation(Dialog(enable = false));
    
protected
  
  final parameter Integer nP = size(coeffs,1) "Number of polynomial coefficients";

equation

  for i in 1:N loop
    
    n_eod[i] = CF_N[i]*sum({coeffs[j]*lambda[i]^(j-1) for j in 1:nP});
  
  end for;

  annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)),
    Documentation(info="<html>
<p>Use of a polynoial based equation for electro-osmotic drag using the reference data as a demonstrable use case.</p>
<h4>References</h4>
<p>Tsai, Shang-Wen, and Yong-Song Chen. &ldquo;A Mathematical Model to Study the Energy Efficiency of a Proton Exchange Membrane Fuel Cell with a Dead-Ended Anode.&rdquo; <i>Applied Energy</i> 188 (February 15, 2017): 151&ndash;59. <a href=\"https://doi.org/10.1016/j.apenergy.2016.11.128\">https://doi.org/10.1016/j.apenergy.2016.11.128</a>. </p>
</html>", revisions="<html>
Copyright &copy; 2004-2024, MODELON AB <br /> The use of this software component is regulated by the licensing conditions for Modelon Libraries. <br />This copyright notice must, unaltered, accompany all components that are derived from, copied from, <br />or by other means have their origin from any Modelon Library.
</html>"));

end Dutta;
