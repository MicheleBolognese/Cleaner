within PEMFCModel.Membrane.Losses.OhmicLoss;

model PEMFCSimplified "PEMFC simplified ohmic loss"
  
  extends PEMFCModel.Membrane.Losses.Templates.BaseOhmicLoss;
  
  parameter .FuelCell.Internal.Units.AreaSpecificResistance R = 1e-6 "Area specific resistance" annotation(Dialog(group = "Input parameters"));

equation
    
  for i in 1:N loop
    
    ASR[i] = R;
    
  end for;

  annotation (Icon(graphics={Text(
          extent={{-64,30},{68,-28}},
          textColor={28,108,200},
          textString="E_ohm")}), Documentation(info="<html>
<p>Ohmic loss expression:</p><p><img src=\"modelica://FuelCell/Resources/images/equations/equation-tdsaN8qv.png\" alt=\"E_loss = -(ASR/A_cell)*i\"/></p>
<p>where parameter ASR represents area specific resistance, A<sub>cell</sub> the membrane area and i the current. </p>
</html>", revisions="<html>
Copyright &copy; 2004-2024, MODELON AB <br /> The use of this software component is regulated by the licensing conditions for Modelon Libraries. <br />This copyright notice must, unaltered, accompany all components that are derived from, copied from, <br />or by other means have their origin from any Modelon Library.
</html>"));

end PEMFCSimplified;
