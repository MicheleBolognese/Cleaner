within PEMFCModel.Membrane.MassTransport.Templates;

partial model PartialWaterPermeation

  extends PEMFCModel.Membrane.MassTransport.Interfaces.PartialWaterFlow(redeclare package Medium_H2O = .Modelon.Media.PreDefined.Liquids.IncompressibleWater, 
  m_flow = {k[i]*d_H2O[i]/mu_H2O[i]*dp_cell[i]/z*A_cell/N*n_cell for i in 1:N});

  .Modelica.Units.SI.Area k[N] "Water Darcy's permeation coefficient";

protected
  
  .Modelica.Units.SI.Pressure p_H2O[N] = 0.5*(p_an_H2O + p_cath_H2O) "Average H2O partial pressure in each membrane node";
  .Modelica.Units.SI.PressureDifference dp_cell[N] = p_an_H2O - p_cath_H2O "H2O partial pressure difference across the cell in each membrane node";
   Medium_H2O.ThermodynamicState state_H2O[N] = Medium_H2O.setState_pTX(p_H2O, T_cell, {1}) "Thermodynamic state of water in each membrane node";
  .Modelica.Units.SI.Density d_H2O[N] = Medium_H2O.density(state_H2O) "Water density";
  .Modelica.Units.SI.DynamicViscosity mu_H2O[N] = Medium_H2O.dynamicViscosity(state_H2O) "Water dynamic viscosity";

  annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={Rectangle(
          extent={{-100,100},{100,-100}},
          lineColor={0,0,0},
          fillPattern=FillPattern.Solid,
          fillColor={255,255,255}), Text(
          extent={{-40,20},{40,-20}},
          textColor={0,0,0},
          textString="WP")}), Diagram(coordinateSystem(preserveAspectRatio=false)),
    Documentation(info="<html>
<p>Template class for the transfer of water across a membrane due to pressure differences on either side of the membrane (i.e., permeation).</p>
</html>", revisions="<html>
Copyright &copy; 2004-2024, MODELON AB <br /> The use of this software component is regulated by the licensing conditions for Modelon Libraries. <br />This copyright notice must, unaltered, accompany all components that are derived from, copied from, <br />or by other means have their origin from any Modelon Library.
</html>"));

end PartialWaterPermeation;
