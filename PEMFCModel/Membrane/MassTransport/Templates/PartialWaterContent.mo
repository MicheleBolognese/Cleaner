within PEMFCModel.Membrane.MassTransport.Templates;

partial model PartialWaterContent

  extends PEMFCModel.Membrane.MassTransport.Interfaces.PartialWaterModel(enableInternal = false);

  final parameter .Modelica.Units.SI.MolarVolume v_dry = EW_m/rho_dry_m "Dry membrane molar-based specific volume" annotation(Dialog(enable = false));
  constant .Modelica.Units.SI.MolarVolume v_w = 1.81e-5 "Water molar-based specific volume" annotation(Dialog(enable = false));
  
  .Modelica.Units.SI.DimensionlessRatio f_w[N] "Volume fraction of water in each membrane node";
  Real lambda[N](each min = 0) "Membrane average water content in each node | Number of water molecules per available site (HSO3- group)";
  Real lambda_an[N](each min = 0) "Membrane water content at the anode-side in each node";
  Real lambda_cath[N](each min = 0) "Membrane water content at the cathode-side in each node";
  
equation

    for i in 1:N loop
      
      f_w[i] = lambda[i]*v_w/(v_dry + lambda[i]*v_w);
    
    end for;
  
  annotation (defaultComponentName="membraneWaterProps",Icon(coordinateSystem(preserveAspectRatio=false), graphics={Rectangle(
          extent={{-100,100},{100,-100}},
          lineColor={0,0,0},
          fillPattern=FillPattern.Solid,
          fillColor={255,255,255}), Text(
          extent={{-40,20},{40,-20}},
          textColor={0,0,0},
          textString="WC")}),                                                                              Diagram(coordinateSystem(preserveAspectRatio=false)),
    Documentation(revisions="<html>
Copyright &copy; 2004-2024, MODELON AB <br /> The use of this software component is regulated by the licensing conditions for Modelon Libraries. <br />This copyright notice must, unaltered, accompany all components that are derived from, copied from, <br />or by other means have their origin from any Modelon Library.
</html>", info="<html>
<p>Template class for the determination of the water content of a membrane. </p>
</html>"));

end PartialWaterContent;
