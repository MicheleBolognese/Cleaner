within PEMFCModel.Membrane.MassTransport.Templates;

partial model PartialElectroOsmoticDrag

  extends PEMFCModel.Membrane.MassTransport.Interfaces.PartialWaterFlow( m_flow = {if p_an_H2O[i] < .Modelica.Constants.eps then 0 else n_eod[i]*I_cell[i]/.Modelica.Constants.F*MM_H2O*n_cell for i in 1:N});

   input .Modelica.Units.SI.Current I_cell[N] "Electrical current per each cell node" annotation(Dialog(tab = "Internal", enable = enableInternal, group = "Input variables"));
  .Modelica.Units.SI.DimensionlessRatio n_eod[N] "Electro-osmotic drag coefficient (approximately 0.27 to 5)";

  annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={Rectangle(
          extent={{-100,100},{100,-100}},
          lineColor={0,0,0},
          fillPattern=FillPattern.Solid,
          fillColor={255,255,255}), Text(
          extent={{-40,20},{40,-20}},
          textColor={0,0,0},
          textString="EOD")}), Diagram(coordinateSystem(preserveAspectRatio=false)),
    Documentation(info="<html>
<p>Template class for defining the transfer of water across a membrane due to proton exchange across the membrane.</p>
</html>", revisions="<html>
Copyright &copy; 2004-2024, MODELON AB <br /> The use of this software component is regulated by the licensing conditions for Modelon Libraries. <br />This copyright notice must, unaltered, accompany all components that are derived from, copied from, <br />or by other means have their origin from any Modelon Library.
</html>"));

end PartialElectroOsmoticDrag;
