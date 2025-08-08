within PEMFCModel.Membrane.MassTransport.GasDiffusion;

model ArrheniusPolyWaterFraction "Arrhenius: K_perm = exp(Ea/R*(1/Tref-1/T) with A = f(lambda) where f is a polynomial function"
    
  extends PEMFCModel.Membrane.MassTransport.Templates.PartialGasDiffusion;

  parameter Boolean useT_ref = true "= false to remove 1/T_ref from Arrhenius equation" annotation(Dialog(group = "Input parameters"));

  parameter .Modelica.Units.SI.MolarEnergy E[nDiff] = fill(19830, nDiff) "Kind of activation energy for each chemical species" annotation(Dialog(group = "Input parameters"));
  parameter .Modelica.Units.SI.Temperature T_ref[nDiff] = fill(303, nDiff) "Reference temperature for each chemiacl species" annotation (Dialog(enable = useT_ref, group = "Input parameters"));

  parameter Real coeffs[nDiff,:] = fill(0, nDiff, 1) "Coefficients for each chemical species: c[1] + c[2]*f_w + ..." annotation(Dialog(group = "Input parameters"));

protected
  
  final parameter Integer nP = size(coeffs, 2) "Number of polynomial coefficients";
  Real A[N, nDiff] "Arrhenius pre-exponential factor for each chemical species in each node";

equation

  for i in 1:N loop
    
    for j in 1:nDiff loop
    
      A[i, j] = sum({coeffs[j, k]*f_w[i]^(k - 1) for k in 1:nP});
    
    end for;

  end for;
  if useT_ref then

    for i in 1:N loop

      for j in 1:nDiff loop

        K_perm[i, j] = CF_N[i, j]*A[i, j]*exp(E[j]/.Modelica.Constants.R*(1/T_ref[j] - 1/T_cell[i]));

      end for;

    end for;

  else

    for i in 1:N loop

      for j in 1:nDiff loop

        K_perm[i, j] = CF_N[i, j]*A[i, j]*exp(-E[j]/(.Modelica.Constants.R*T_cell[i]));

      end for;

    end for;

  end if;

  annotation (
    Icon(coordinateSystem(preserveAspectRatio=false)),
    Diagram(coordinateSystem(preserveAspectRatio=false)),
    Documentation(info="<html>
<p>A model for applying the Arrhenius equation using a polynomial volume fraction of water in the membrane to calculate the pre-exponential factor.</p>
<p>To use:</p>
<ol>
<li>Specify and array of strings of the gas species of interest. The species MUST be present in the media on both sides of the membrane.</li>
<li>Provide array of exponential values for each species.</li>
<li>Provide array of reference temperature values.</li>
<li>Provide a matrix of coefficients for each species. Note that any number of coefficients are accepted. If all species do not have the specified coefficient, zero can be used.</li>
</ol>
<p><u>Examples using data from Tsai et. al.:</u></p>
<ul>
<li>Example for N2: </li>
<p>diffusiveSpecies={&quot;N2&quot;}</p>
<p>E ={19830}</p>
<p>T_ref={303}</p>
<p>coeffs = {{0.0295,1.21,-1.93}*1e-14}</p>
<li>Example for H2, O2, and N2: </li>
<p>diffusiveSpecies={&quot;H2&quot;,&quot;O2&quot;,&quot;N2&quot;}</p>
<p>E ={18000,20000,19830}</p>
<p>T_ref=fill(303,&nbsp;nDiff)</p>
<p>coeffs = {{0.29,2.2,0}*1e-14,{0.11,1.9,0}*1e-14,{0.0295,1.21,-1.93}*1e-14};</p>
</ul>
<h4>References</h4>
<p>Tsai, Shang-Wen, and Yong-Song Chen. &ldquo;A Mathematical Model to Study the Energy Efficiency of a Proton Exchange Membrane Fuel Cell with a Dead-Ended Anode.&rdquo; <i>Applied Energy</i> 188 (February 15, 2017): 151&ndash;59. <a href=\"https://doi.org/10.1016/j.apenergy.2016.11.128\">https://doi.org/10.1016/j.apenergy.2016.11.128</a>. </p>
</html>", revisions="<html>
Copyright &copy; 2004-2024, MODELON AB <br /> The use of this software component is regulated by the licensing conditions for Modelon Libraries. <br />This copyright notice must, unaltered, accompany all components that are derived from, copied from, <br />or by other means have their origin from any Modelon Library.
</html>"));

end ArrheniusPolyWaterFraction;
