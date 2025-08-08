within PEMFCModel.Membrane.Losses.ActivationLoss;

model PEMFCElectroChem "PEMFC electrochemical activation loss component "
  
  extends PEMFCModel.Membrane.Losses.Templates.BaseActivationLoss;

  parameter Real alpha = 0.5 "Charge transfer coefficient" annotation(Dialog(group = "Input parameters"));
  parameter Real n_e = 2 "Number of exchanged electrons" annotation(Dialog(group = "Input parameters"));
  parameter .Modelica.Units.SI.Thickness z = 20e-6 "Membrane thickness" annotation(Dialog(enable = enableInternal, tab = "Internal", group = "Input parameters"));
  parameter .Modelica.Units.SI.CurrentDensity j0_ref = 5e-5 "Reference exchange current density for activation loss" annotation(Dialog(group = "Input parameters"));
  parameter Real Rc = 100 "Roughness factor of the cathode, i.e. catalyst active area per unit cathode area" annotation(Dialog(group = "Input parameters"));
  parameter Real gamma_ORR = 0.5 "O2 reduction reaction dependency factor on O2 partial pressure" annotation(Dialog(group = "Input parameters"));
  parameter .FuelCell.Internal.Units.ActivationEnergy Ea_ORR = 66e3 "Activation energy for ORR on the catalyst" annotation(Dialog(group = "Input parameters"));
  parameter .FuelCell.Internal.Units.ActivationEnergy Ea_k_H2 = 21.03e3 "Activation energy for H2 permeation" annotation(Dialog(group = "Input parameters"));
  parameter Real A_k_H2 = 6.6e-11 "Pre-exponetial factor for H2 permeation" annotation(Dialog(group = "Input parameters"));
  constant .Modelica.Units.SI.Temperature Tref_25C = 298.15 "Reference temperature at 25°C" annotation(Dialog(enable = false));

  input .Modelica.Units.SI.Temperature T_cell[N] "Cell substrate temperature in each node" annotation(Dialog(enable = enableInternal, tab = "Internal", group = "Input variables"));
  input .Modelica.Units.SI.Pressure pH2[N] "H2 partial pressure on membrane anode surface in each node" annotation(Dialog(enable = enableInternal, tab = "Internal", group = "Input variables"));
  input .Modelica.Units.SI.Pressure pO2[N] "O2 partial pressure on membrane cathode surface in each node" annotation(Dialog(enable = enableInternal, tab = "Internal", group = "Input variables"));

  .Modelica.Units.SI.CurrentDensity j_0[N] "Exchange current density in each node";
  .FuelCell.Internal.Units.ReactantPermeability k_H2[N] "H2 permeability in each node";
  .Modelica.Units.SI.CurrentDensity j_loss[N](each min = 0) "Current density loss in each node";

equation
 
  for i in 1:N loop
    
    j_0[i] = j0_ref*Rc*(pO2[i]/.FuelCell.Internal.Units.pNTP)^gamma_ORR*exp(-Ea_ORR/(.Modelica.Constants.R*T_cell[i])*(1-T_cell[i]/Tref_25C));
    k_H2[i] = A_k_H2*exp(-Ea_k_H2/(.Modelica.Constants.R*T_cell[i]));
    j_loss[i] = n_e*k_H2[i]*pH2[i]*.FuelCell.Internal.Units.F/z;
    E_loss_cell[i] = (.Modelica.Constants.R*T_cell[i]/(n_e*alpha*.FuelCell.Internal.Units.F)*.Modelica.Math.log((j_ionic[i]+j_loss[i])/j_0[i]));
    
  end for;

  annotation (Icon(graphics={Text(
          extent={{-64,30},{68,-28}},
          textColor={28,108,200},
          textString="E_act")}), Documentation(info="<html>
<p>Activation loss is parameterized by number of exchanged electrons <span style=\"font-family: Courier New;\">n_e</span> and the charge transfer coefficient <span style=\"font-family: Courier New;\">alpha</span>.</p>
<p><img src=\"modelica://FuelCell/Resources/images/equations/equation-EIAi0oHN.png\" alt=\"E_loss = RT_cell/(alpha*n_e*F)*ln((j_Ionic+j_n)/j_0)\"/></p>
<p>The required inputs are the substrate temperature T_cell, pressure of H2 and O2 on anode and cathode surfaces respectively, exchange current density j_0, current density loss j_loss and the effective current j_Ionic.</p>
<p>The current density loss and exchange current density are calculated as follows:</p>
<p><img src=\"modelica://FuelCell/Resources/images/equations/equation-Fr8eUaVR.png\" alt=\"j_loss = n_e*k_H2*pH2*F/z\"/></p>
<p><br>where: </p>
<p><img src=\"modelica://FuelCell/Resources/images/equations/equation-LDZyoeRy.png\" alt=\"k_H2 = Ak_H2*exp(-Eak_H2/(R*T_cell))\"/></p>
<p><br>and</p>
<p><img src=\"modelica://FuelCell/Resources/images/equations/equation-jU9K8Iyl.png\" alt=\"j_0 = j_0ref*R_c*(pO2/p_atm)^gamma_ORR * exp(-Ea_ORR/(R*T_cell)*(1-T_cell/T_ref25C))\"/></p>
</html>", revisions="<html>
Copyright &copy; 2004-2024, MODELON AB <br /> The use of this software component is regulated by the licensing conditions for Modelon Libraries. <br />This copyright notice must, unaltered, accompany all components that are derived from, copied from, <br />or by other means have their origin from any Modelon Library.
</html>"));

end PEMFCElectroChem;
