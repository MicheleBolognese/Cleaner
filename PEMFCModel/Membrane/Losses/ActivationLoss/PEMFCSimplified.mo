within PEMFCModel.Membrane.Losses.ActivationLoss;

model PEMFCSimplified "PEMFC simplified activation loss"
  
  extends PEMFCModel.Membrane.Losses.Templates.BaseActivationLoss;

  parameter Real a = 0.0302 "Tafel slope" annotation(Dialog(group = "Input parameters"));
  parameter .Modelica.Units.SI.CurrentDensity j0 = 1 "Exchange current density" annotation(Dialog(group = "Input parameters"));

equation
 
  for i in 1:N loop
    
    E_loss_cell[i] = a*.Modelica.Math.log(max(j_ionic[i]/j0, .Modelica.Constants.eps) + 1);
 
  end for;

  annotation (Icon(graphics={Text(
          extent={{-64,30},{68,-28}},
          textColor={28,108,200},
          textString="E_act")}), Documentation(info="<html>
<p>Activation loss based on the Butler-Volmer expression simplified to the Tafel expression:</p>
<p><img src=\"modelica://FuelCell/Resources/images/equations/equation-9TONjlsH.png\" alt=\"E_loss = a * ln (j/j0 +1)\"/></p>
<p>where <i>a</i> represents the Tafel slope, <i>j</i> the current density and <i>j0</i> the so-called exchange current density. </p>
</html>", revisions="<html>
Copyright &copy; 2004-2024, MODELON AB <br /> The use of this software component is regulated by the licensing conditions for Modelon Libraries. <br />This copyright notice must, unaltered, accompany all components that are derived from, copied from, <br />or by other means have their origin from any Modelon Library.
</html>"));

end PEMFCSimplified;
