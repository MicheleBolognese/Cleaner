within PEMFCModel.Membrane.Losses.ActivationLoss;
model PEMFCEmpirical "PEMFC empirical activation loss"
  
  extends PEMFCModel.Membrane.Losses.Templates.BaseActivationLoss;

  parameter Real alpha = 0.5 "Charge transfer coefficient" annotation(Dialog(group = "Input parameters"));
  constant Real n_e = 2.0 "Number of exchanged electrons" annotation(Dialog(enable = false));
  parameter .Modelica.Units.SI.CurrentDensity j_loss = 100 "Current density loss" annotation(Dialog(group = "Input parameters"));
  parameter .Modelica.Units.SI.CurrentDensity j_0 = 1 "Exchange current density" annotation(Dialog(group = "Input parameters"));
  input .Modelica.Units.SI.Temperature T_cell[N] "Cell substrate temperature in each node" annotation(Dialog(enable = enbaleInternal, tab = "Internal", group = "Input variables"));

equation
 
  for i in 1:N loop
    
    E_loss_cell[i] = .Modelica.Constants.R*T_cell[i]/(n_e*alpha*.FuelCell.Internal.Units.F)*Modelica.Math.log(max((j_ionic[i] + j_loss)/j_0, Modelica.Constants.eps)); //Modelon.Math.Smoothing.above(.Modelica.Constants.R*T_cell[i]/(n_e*alpha*.FuelCell.Internal.Units.F)*Modelica.Math.log(max((j_ionic[i] + j_loss)/j_0, Modelica.Constants.eps)),0,0.5);
 
  end for;
  
  annotation (Icon(graphics={Text(
          extent={{-64,30},{68,-28}},
          textColor={28,108,200},
          textString="E_act")}), Documentation(info="<html>
<p>Activation loss is parameterized by number of exchanged electrons <span style=\"font-family: Courier New;\">n_e</span> and the charge transfer coefficient <span style=\"font-family: Courier New;\">alpha</span>.</p>
<p><img src=\"modelica://FuelCell/Resources/images/equations/equation-EIAi0oHN.png\" alt=\"E_loss = RT_cell/(alpha*n_e*F)*ln((j_Ionic+j_n)/j_0)\"/></p>
<p>The required inputs are the substrate temperature T_cell, exchange current density j_0, current density loss j_n and the effective current j_ionic.</p>
</html>", revisions="<html>
Copyright &copy; 2004-2024, MODELON AB <br /> The use of this software component is regulated by the licensing conditions for Modelon Libraries. <br />This copyright notice must, unaltered, accompany all components that are derived from, copied from, <br />or by other means have their origin from any Modelon Library.
</html>"));

end PEMFCEmpirical;