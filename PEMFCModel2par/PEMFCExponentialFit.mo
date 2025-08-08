within PEMFCModel2par;

model PEMFCExponentialFit "PEMFC exponentially-fitted concentration loss"
  
  extends .PEMFCModel.Membrane.Losses.Templates.BaseConcentrationLoss;

  input .Modelica.Units.SI.Voltage m (start = 3e-4) "Pre-exponential coefficient" annotation(Dialog(group = "Input parameters"));
  input Real n(unit = "m2/A", start = 3.2e-4) "Exponential coefficient" annotation(Dialog(group = "Input parameters"));

equation
    
 for i in 1:N loop
    
  E_loss_cell[i] = m*exp(n*j_ionic[i]);
    
 end for;

 annotation (Icon(graphics={Text(
          extent={{-64,34},{68,-24}},
          textColor={28,108,200},
          textString="E_conc")}), Documentation(info="<html>
<p>Concentration loss is parameterized by the experimentally derived values of m and n coefficients.</p>
<p><img src=\"modelica://FuelCell/Resources/images/equations/equation-NsFpmSpW.png\" alt=\"E_loss = m*exp(n*j)\"/></p>
</html>", revisions="<html>
Copyright &copy; 2004-2024, MODELON AB <br /> The use of this software component is regulated by the licensing conditions for Modelon Libraries. <br />This copyright notice must, unaltered, accompany all components that are derived from, copied from, <br />or by other means have their origin from any Modelon Library.
</html>"));

end PEMFCExponentialFit;
