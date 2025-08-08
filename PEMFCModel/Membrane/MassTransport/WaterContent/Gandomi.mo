within PEMFCModel.Membrane.MassTransport.WaterContent;
model Gandomi "Correlation by Gandomi for membrane water content"
 
  extends PEMFCModel.Membrane.MassTransport.Templates.PartialWaterContent;
  
//protected

  Real a_an_H2O[N] = p_an_H2O./Medium_H2O.saturationPressure_TX(T_an, {1}) "H2O activity at anode side in each node";
  Real a_cath_H2O[N] = p_cath_H2O./Medium_H2O.saturationPressure_TX(T_cath, {1}) "H2O activity at cathode side in each node";
  Real a_cell_H2O[N] = 0.5*(a_an_H2O + a_cath_H2O) "Membrane average water activity in each node";

equation
  
  for i in 1:N loop
      
    lambda[i] = CF_N[i]*.PEMFCModel.Membrane.MassTransport.WaterContent.WaterContentGandomiFunction(a = a_cell_H2O[i]);
    lambda_an[i] = CF_N[i]*.PEMFCModel.Membrane.MassTransport.WaterContent.WaterContentGandomiFunction(a = a_an_H2O[i]);
    lambda_cath[i] = CF_N[i]*.PEMFCModel.Membrane.MassTransport.WaterContent.WaterContentGandomiFunction(a = a_cath_H2O[i]);
  
  end for;

  annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)),
    Documentation(info="<html>
<p>Membrane water content as a function of the activity of water as defined by Gandomi et. al.</p>
<h4>References</h4>
<p>Gandomi, Yasser Ashraf, M. D. Edmundson, F. C. Busby, and Matthew M. Mench. &ldquo;Water Management in Polymer Electrolyte Fuel Cells through Asymmetric Thermal and Mass Transport Engineering of the Micro-Porous Layers.&rdquo; <i>Journal of The Electrochemical Society</i> 163, no. 8 (June 16, 2016): F933. <a href=\"https://doi.org/10.1149/2.1331608jes\">https://doi.org/10.1149/2.1331608jes</a>. </p>
</html>", revisions="<html>
Copyright &copy; 2004-2024, MODELON AB <br /> The use of this software component is regulated by the licensing conditions for Modelon Libraries. <br />This copyright notice must, unaltered, accompany all components that are derived from, copied from, <br />or by other means have their origin from any Modelon Library.
</html>"));

end Gandomi;