within PEMFCModel.Membrane.Losses.OhmicLoss;
model PEMFCElectroChem "PEMFC electrochemical ohmic loss"
    
  extends PEMFCModel.Membrane.Losses.Templates.BaseOhmicLoss;

  parameter .Modelica.Units.SI.Thickness z = 20e-6  "Membrane thickness" annotation(Dialog(enable = enableInternal, tab = "Internal", group = "Input parameters"));
  input .Modelica.Units.SI.Temperature T_cell[N] "Cell substrate temperature in each node" annotation(Dialog(enable = enableInternal, tab = "Internal", group = "Input variables"));
  input Real lambda_mean[N] "Average water content in each membrane node" annotation(Dialog(enable = enableInternal, tab = "Internal", group = "Input variables"));
    
  .Modelica.Units.SI.Conductivity sigma[N] "Membrane conductivity in each node";

equation
    
  for i in 1:N loop
    
    sigma[i] = (5.139e-3*lambda_mean[i]-3.26e-3)*exp(1268*(1/303-1/T_cell[i]))*1e2;
    ASR[i] = z/max(.Modelica.Constants.eps, sigma[i]);
  
  end for;

  annotation (Icon(graphics={Text(
          extent={{-64,30},{68,-28}},
          textColor={28,108,200},
          textString="E_ohm")}), Documentation(info="<html>
<p>Ohmic loss expression:</p><p><img src=\"modelica://FuelCell/Resources/images/equations/equation-tdsaN8qv.png\" alt=\"E_loss = -(ASR/A_cell)*i\"/></p>
<p>where ASR represents area specific resistance, A<sub>cell</sub> the membrane area and i the current. </p>
<p>Humidity effects at both the anode and cathode are included in the calculation of internal resistance (ASR) given as:</p>
<p><img src=\"modelica://FuelCell/Resources/images/equations/equation-AehaYlYY.png\" alt=\"ASR = z/sigma\"/></p>
<p>where z is the membrane thickness and &sigma; is the membrane conductivity, based on the water content &lambda; and cell temperature, fitted as follows:</p>
<p><img src=\"modelica://FuelCell/Resources/images/equations/equation-JZ9tzsiN.png\" alt=\"sigma = (0.5139*lambda-0.326)*exp(1268*(1/303 - 1/T))\"/></p>
<p>Pressure of H2O on anode and cathode surface of membrane are required as input.</p>
<p>A simple average of water vapor activities at the anode and cathode is then used to calculate the average water content of the membrane.</p>
</html>", revisions="<html>
Copyright &copy; 2004-2024, MODELON AB <br /> The use of this software component is regulated by the licensing conditions for Modelon Libraries. <br />This copyright notice must, unaltered, accompany all components that are derived from, copied from, <br />or by other means have their origin from any Modelon Library.
</html>"));

end PEMFCElectroChem;