within PEMFCModel.Membrane.Losses.ConcentrationLoss;

model PEMFCElectroChem "PEMFC electrochemical concentration loss"
  
  extends PEMFCModel.Membrane.Losses.Templates.BaseConcentrationLoss;

  parameter Real alpha = 0.5 "Charge transfer coefficient" annotation(Dialog(group = "Input parameters"));
  parameter Real n_e = 2 "Number of exchanged electrons" annotation(Dialog(group = "Input parameters"));
  parameter .Modelica.Units.SI.CurrentDensity j_lim = 20000 "Limiting exchange current density" annotation(Dialog(group = "Input parameters"));
  
  input .Modelica.Units.SI.Temperature[N] T_cell "Cell substrate temperature in each node" annotation(Dialog(enable = enableInternal, tab = "Internal", group = "Input variables"));

equation
    
  for i in 1:N loop
    
    E_loss_cell[i] = -.Modelica.Constants.R*T_cell[i]/(n_e*.FuelCell.Internal.Units.F)*(1+1/alpha)*.Modelica.Math.log(max(.Modelica.Constants.eps, 1.0 - j_ionic[i]/j_lim));
    
  end for;

  annotation (Icon(graphics={Text(
          extent={{-64,34},{68,-24}},
          textColor={28,108,200},
          textString="E_conc")}), Documentation(info="<html>
<p>Concentration loss is parameterized by number of exchanged electrons <span style=\"font-family: Courier New;\">n_e</span> and the charge transfer coefficient <span style=\"font-family: Courier New;\">alpha</span>.</p>
<p><img src=\"modelica://FuelCell/Resources/images/equations/equation-KoYyuB5v.png\" alt=\"E_loss = (RT_cell/(n_e*F))*(1+ 1/alpha)*ln(j_L/(j_L-j_Ionic))\"/></p>
<p>The required inputs are the substrate temperature T_cell, maximum exchange current density j_L and the effective current j_ionic.</p>
</html>", revisions="<html>
Copyright &copy; 2004-2024, MODELON AB <br /> The use of this software component is regulated by the licensing conditions for Modelon Libraries. <br />This copyright notice must, unaltered, accompany all components that are derived from, copied from, <br />or by other means have their origin from any Modelon Library.
</html>"));

end PEMFCElectroChem;
