within PEMFCModel.Membrane.MassTransport.Templates;

partial model PartialGasDiffusion

  extends PEMFCModel.Membrane.MassTransport.Interfaces.PartialTransport;

  parameter Real CF_0[nDiff] = ones(nDiff) "Correction factor (1 by default)" annotation (Dialog(tab = "Advanced", group = "Input parameters"));
  parameter Real CF_N[N, nDiff] = {{CF_0[j] for j in 1:nDiff} for i in 1:N} "Set if CF_0 is non-uniform" annotation(Dialog(tab = "Advanced", group = "Input parameters"));

  input Real lambda[N] "Membrane average water content in each node" annotation(Dialog(tab = "Internal", enable = enableInternal, group = "Input variables"));
  input .Modelica.Units.SI.DimensionlessRatio f_w[N] "Volume fraction of water in each membrane node" annotation (Dialog(tab = "Internal", enable = enableInternal, group = "Input variables"));

  parameter String diffusiveSpecies[:] = fill("", 0) "Name of chemical species to include in diffusion (e.g., {'N2'})" annotation(Dialog(group = "Input parameters"));
  final parameter Integer nDiff = size(diffusiveSpecies, 1);

  Real K_perm[N, nDiff](each unit = "mol/(s.m.Pa)") "Permeability coefficient (diffusion*solubility)";

  .Modelica.Units.SI.MassFlowRate mX_flow_diff[N, nDiff] = {{-K_perm[i, j]*(p_cath_partial[i, iDiffusive_cath[j]] - p_an_partial[i, iDiffusive_an[j]])/z*Medium_an.MMX[iDiffusive_an[j]]*A_cell/N*n_cell for j in 1:nDiff} for i in 1:N} "Mass flow rates of diffusive species only (> 0 from anode to cathode)";

  .Modelica.Units.SI.MassFlowRate mX_flow_an[N, Medium_an.nS] "Diffusive mass flow rates of each chemical species in each node (> 0 from anode to cathode)";
  .Modelica.Units.SI.MassFlowRate mX_flow_cath[N, Medium_cath.nS] "Diffusive mass flow rates of each chemical species in each node (> 0 from cathode to anode)";

protected
  
  final parameter Integer iDiffusive_an[nDiff] = Medium_an.substanceIndexVector(diffusiveSpecies) "Indices of diffusive gas species in the anodic medium";
  final parameter Integer iDiffusive_cath[nDiff] = Medium_cath.substanceIndexVector(diffusiveSpecies) "Indices of diffusive gas species in the cathodic medium";
  final parameter Integer iMapDiff_an[Medium_an.nS] = {.Modelica.Math.Vectors.find(i, iDiffusive_an) for i in 1:Medium_an.nS};
  final parameter Integer iMapDiff_cath[Medium_cath.nS] = {.Modelica.Math.Vectors.find(i, iDiffusive_cath) for i in 1:Medium_cath.nS};

equation

  for i in 1:N loop
    
    for j in 1:Medium_an.nS loop
    
      mX_flow_an[i, j] = if iMapDiff_an[j] == 0 then 0 else mX_flow_diff[i, iMapDiff_an[j]];
    
    end for;
    for j in 1:Medium_cath.nS loop
      
      mX_flow_cath[i, j] = if iMapDiff_cath[j] == 0 then 0 else -mX_flow_diff[i, iMapDiff_cath[j]];
    
    end for;
  
  end for;

  annotation (Icon(coordinateSystem(preserveAspectRatio=false), graphics={Rectangle(
          extent={{-100,100},{100,-100}},
          lineColor={0,0,0},
          fillPattern=FillPattern.Solid,
          fillColor={255,255,255}), Text(
          extent={{-40,20},{40,-20}},
          textColor={0,0,0},
          textString="GD")}), Diagram(coordinateSystem(preserveAspectRatio=false)),
    Documentation(info="<html>
<p>Template class for the definition of diffusion of gasses across a membrane. The default equation is based on the partial pressure difference across the membrane (i.e., Fick&apos;s law).</p>
</html>", revisions="<html>
Copyright &copy; 2004-2024, MODELON AB <br /> The use of this software component is regulated by the licensing conditions for Modelon Libraries. <br />This copyright notice must, unaltered, accompany all components that are derived from, copied from, <br />or by other means have their origin from any Modelon Library.
</html>"));

end PartialGasDiffusion;
