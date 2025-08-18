within SAMU;

partial model PartialCellTransport_1 "General fuel cell membrane interface with additional transport model"

  //extends .PEMFCModel.Membrane.Templates.PartialCell;

  parameter .Modelica.Units.SI.Length z = 20e-6 "Membrane thickness" annotation (Dialog(enable = enable_setting, group = "Membrane characteristics"));
  parameter .Modelica.Units.SI.Density rho_dry_m = 1980 "Dry membrane density" annotation (Dialog(enable = enable_setting, group = "Membrane characteristics"));
  parameter .Modelica.Units.SI.MolarMass EW_m = 1.1 "Membrane equivalent weight" annotation (Dialog(enable = enable_setting, group = "Membrane characteristics"));
  
  parameter String diffusiveSpecies[:] = fill("",0) "Species diffusing through membrane" annotation(Dialog(enable = enable_setting, group = "Membrane characteristics"));

  .Modelica.Units.SI.Pressure p_an_partial[N, nS_an] = {{max(.Modelica.Constants.eps, p_an[i]*y_an[i, j]) for j in 1:nS_an} for i in 1:N} "Partial pressure of each chemical species in each node at membrane-anode interface ";
  .Modelica.Units.SI.Pressure p_cath_partial[N, nS_cath] = {{max(.Modelica.Constants.eps, p_cath[i]*y_cath[i, j]) for j in 1:nS_cath} for i in 1:N} "Partial pressure of each chemical species in each node at membrane-cathode interface";

  replaceable model WaterContent = .PEMFCModel.Membrane.MassTransport.WaterContent.Gandomi constrainedby .PEMFCModel.Membrane.MassTransport.Templates.PartialWaterContent(final enableInternal = false) annotation (Dialog(tab = "Mass transport", group = "General properties",enable = true));
  
  parameter Real CF_0_waterContent = 1.0 "Calibration factor for water content" annotation(Dialog(tab = "Mass transport", group = "General properties",enable =  true));
  parameter Real CF_N_waterContent[N] = fill(CF_0_waterContent, N) "Calibration factor for water content (set if non-uniform)" annotation(Dialog(tab = "Mass transport", group = "General properties",enable = true));
  
  WaterContent waterContent(
    redeclare package Medium_an = Medium_an,
    redeclare package Medium_cath = Medium_cath,
    final N = N,
    final n_cell = n_cell,
    final A_cell = A_cell,
    final z = z,
    final rho_dry_m = rho_dry_m,
    final EW_m = EW_m,
    final p_an_partial = p_an_partial,
    final p_cath_partial = p_cath_partial,
    final T_an = T_an,
    final T_cath = T_cath,
    final T_cell = T_cell,
    final iH2O_an = iH2O_an[1],
    final iH2O_cath = iH2O_cath[1],
    final CF_0 = CF_0_waterContent,
    final CF_N = CF_N_waterContent) annotation (Placement(transformation(extent={{-108.75609436022503,94.04853189230951},{-88.75609436022503,114.04853189230951}},rotation = 0.0,origin = {0.0,0.0})));
  
  replaceable model WaterDiffusion = .PEMFCModel.Membrane.MassTransport.WaterDiffusion.ZeroFlow constrainedby .PEMFCModel.Membrane.MassTransport.Templates.PartialWaterDiffusion(final enableInternal = false) "Water diffusion through membrane" annotation(choicesAllMatching = true, Dialog(tab = "Mass transport", group = "Water transport",enable = true));
  
  parameter Real CF_0_waterDiffusion = 1.0 "Calibration factor for water diffusion" annotation(Dialog(tab = "Mass transport", group = "Water transport",enable =  true));
  parameter Real CF_N_waterDiffusion[N] = fill(CF_0_waterDiffusion, N) "Calibration factor for water diffusion (set if non-uniform)" annotation(Dialog(tab = "Mass transport", group = "Water transport",enable = true));  

  WaterDiffusion waterDiffusion(
    redeclare package Medium_an = Medium_an,
    redeclare package Medium_cath = Medium_cath,
    final N = N,
    final n_cell = n_cell,
    final A_cell = A_cell,
    final z = z,
    final rho_dry_m = rho_dry_m,
    final EW_m = EW_m,
    final p_an_partial = p_an_partial,
    final p_cath_partial = p_cath_partial,
    final T_an = T_an,
    final T_cath = T_cath,
    final T_cell = T_cell,
    final iH2O_an = iH2O_an[1],
    final iH2O_cath = iH2O_cath[1],
    final lambda = waterContent.lambda,
    final lambda_an = waterContent.lambda_an,
    final lambda_cath = waterContent.lambda_cath,
    final CF_0 = CF_0_waterDiffusion,
    final CF_N = CF_N_waterDiffusion) annotation (Placement(transformation(extent={{-88.7413087329926,94.04853189230951},{-68.7413087329926,114.04853189230951}},rotation = 0.0,origin = {0.0,0.0})));

  replaceable model ElectroOsmoticDrag = .PEMFCModel.Membrane.MassTransport.ElectroOsmoticDrag.ZeroFlow constrainedby .PEMFCModel.Membrane.MassTransport.Templates.PartialElectroOsmoticDrag(final enableInternal = false) "Electro-osmotic drag through membrane" annotation (choicesAllMatching = true, Dialog(tab = "Mass transport", group = "Water transport",enable = true));

  parameter Real CF_0_eod = 1.0 "Calibration factor for electro-osmotic drag" annotation(Dialog(tab = "Mass transport", group = "Water transport",enable = true));
  parameter Real CF_N_eod[N] = fill(CF_0_eod, N) "Calibration factor for electro-osmotic drag (set if non-uniform)" annotation(Dialog(tab = "Mass transport", group = "Water transport",enable = true));    

  ElectroOsmoticDrag electroOsmoticDrag(
    redeclare package Medium_an = Medium_an,
    redeclare package Medium_cath = Medium_cath,
    final N = N,
    final n_cell = n_cell,
    final A_cell = A_cell,
    final z = z,
    final rho_dry_m = rho_dry_m,
    final EW_m = EW_m,
    final I_cell =  j_ionic*A_cell/N /*fill(I_stack/N, N) */,
    final p_an_partial = p_an_partial,
    final p_cath_partial = p_cath_partial,
    final T_an = T_an,
    final T_cath = T_cath,
    final T_cell = T_cell,
    final iH2O_an = iH2O_an[1],
    final iH2O_cath = iH2O_cath[1],
    final lambda = waterContent.lambda,
    final CF_0 = CF_0_eod,
    final CF_N = CF_N_eod) annotation (Placement(transformation(extent={{-69.0336325774397,94.04853189230954},{-49.033632577439704,114.04853189230954}},rotation = 0.0,origin = {0.0,0.0})));

  replaceable model GasDiffusion = .PEMFCModel.Membrane.MassTransport.GasDiffusion.ZeroFlow constrainedby .PEMFCModel.Membrane.MassTransport.Templates.PartialGasDiffusion (final enableInternal = false) "Gas permeation through membrane" annotation (choicesAllMatching = true, Dialog(tab = "Mass transport", group = "Gas diffusion",enable = true));
  
  parameter Real CF_0_gasDiffusion[size(diffusiveSpecies,1)] = fill(1.0,size(diffusiveSpecies,1)) "Calibration factor for each diffusive sepcies" annotation(Dialog(tab = "Mass transport", group = "Gas diffusion",enable = true));
  parameter Real CF_N_gasDiffusion[N, size(diffusiveSpecies,1)] = fill(CF_0_gasDiffusion, N) "Calibration factor for each diffusive species (set if non-uniform)" annotation(Dialog(tab = "Mass transport", group = "Gas diffusion",enable = true));  

  GasDiffusion gasDiffusion(
    redeclare package Medium_an = Medium_an,
    redeclare package Medium_cath = Medium_cath,
    final N = N,
    final n_cell = n_cell,
    final A_cell = A_cell,
    final z = z,
    final rho_dry_m = rho_dry_m,
    final EW_m = EW_m,
    final p_an_partial = p_an_partial,
    final p_cath_partial = p_cath_partial,
    final T_an = T_an,
    final T_cath = T_cath,
    final T_cell = T_cell,
    final diffusiveSpecies = diffusiveSpecies,
    final lambda = waterContent.lambda,
    final f_w = waterContent.f_w,
    final CF_0 = CF_0_gasDiffusion,
    final CF_N = CF_N_gasDiffusion) annotation (Placement(transformation(extent={{-48.80090599247191,94.02866613914978},{-28.80090599247191,114.02866613914978}},rotation = 0.0,origin = {0.0,0.0})));

/*
  replaceable model WaterPermeation = .PEMFCModel.Membrane.MassTransport.WaterPermeation.ZeroFlow constrainedby
    .PEMFCModel.Membrane.MassTransport.Templates.PartialWaterPermeation(final enableInternal = false) "Water permeation through membrane" annotation(choicesAllMatching = true, Dialog(tab = "Mass transport", group = "Water transport"));

  WaterPermeation waterPermeation(
    redeclare package Medium_an = Medium_an,
    redeclare package Medium_cath = Medium_cath,
    final N = N,
    final n_cell = n_cell,
    final A_cell = A_cell,
    final z = z,
    final rho_dry_m = rho_dry_m,
    final EW_m = EW_m,
    final p_an_partial = p_an_partial,
    final p_cath_partial = p_cath_partial,
    final T_an = T_an,
    final T_cath = T_cath,
    final T_cell = T_cell,
    final iH2O_an = iH2O_an[1],
    final iH2O_cath = iH2O_cath[1],
    final lambda = waterContent.lambda) annotation (Placement(transformation(extent={{-60,80},{-40,100}})));
*/
  
  final .Modelica.Blocks.Interfaces.RealOutput mX_flow_an_transport[N, Medium_an.nS](each unit = "kg/s") "Water/gas transport at anode side (> 0 from anode to cathode)";
  final .Modelica.Blocks.Interfaces.RealOutput mX_flow_cath_transport[N, Medium_cath.nS](each unit = "kg/s") "Water/gas transport at cathode side (> 0 from cathode to anode)";
  final parameter Integer iH2O_an[1] = Medium_an.substanceIndexVector({"H2O"}) "Index of H2O in the anodic medium";
  final parameter Integer iH2O_cath[1] = Medium_cath.substanceIndexVector({"H2O"}) "Index of H2O in the cathodic medium";

protected
  
  final parameter Integer i_required_an[size(diffusiveSpecies, 1)] = Medium_an.substanceIndexVector(diffusiveSpecies) "Indices of required species in the anodic medium";
  final parameter Integer i_required_cath[size(diffusiveSpecies, 1)] = Medium_cath.substanceIndexVector(diffusiveSpecies) "Indices of required species in the cathodic medium";

 initial equation    
  
  for i in 1:size(diffusiveSpecies, 1) loop
  
    assert(i_required_an[i] > 0, "Membrane transport: the selected medium for the anode side does not contain the required substances: see 'required_species'",
    level = AssertionLevel.error);
    assert(i_required_cath[i] > 0, "Membrane transport: the selected medium for the cathode side does not contain the required substances: see 'required_species'",
    level = AssertionLevel.error);
  
  end for;

equation

  // Mass flow rates of each chemical species to mass flow ports per each stack node, excluding reacting species flows
  for i in 1:N loop
    
    for j in 1:nS_an loop
      
      if j == iH2O_an[1] then

        mX_flow_an_transport[i, j] = waterDiffusion.m_flow[i] + electroOsmoticDrag.m_flow[i] /* + waterPermeation.m_flow[i] + gasDiffusion.mX_flow_an[i, j]*/;
      
      else
       
        mX_flow_an_transport[i, j] = gasDiffusion.mX_flow_an[i, j];
      
      end if;
    end for;
    for j in 1:nS_cath loop
      
      if j == iH2O_cath[1] then
        
        mX_flow_cath_transport[i, j] = -waterDiffusion.m_flow[i] - electroOsmoticDrag.m_flow[i] /* - waterPermeation.m_flow[i] + gasDiffusion.mX_flow_cath[i, j] */;
      
      else
      
        mX_flow_cath_transport[i, j] = gasDiffusion.mX_flow_cath[i, j];
      
      end if;
    
    end for;
  
  end for;

  annotation (Documentation(info="<html>
<p>The model includes additional transport properties not associated directly with the ionic consumption/generation in the membrane. More specifically:</p>
<ul>
<li>Water transport due to purely pressure differences across the membrane, concentration gradients, and electro-osmotic drag</li>
<li>Species diffusion across the membrane (e.g., N2)</li>
</ul>
Setting <span style=\"font-family: Courier New;\">useTransport=false</span> will nullify the effect of the water and gaseous diffusions. 
<h4>Parametrization</h4>
<p>The model has parameters determining the electrochemical characteristics, which the user can set. There is also a set of predefined parameters, which the user should not override, as they are common for all membranes of a stack and therefore should be defined for the entire </p>
<h4>References</h4>
<ol>
<li>Rabbani, Abid, and Masoud Rokni. &ldquo;Effect of Nitrogen Crossover on Purging Strategy in PEM Fuel Cell Systems.&rdquo; <i>Applied Energy</i> 111 (November 1, 2013): 1061&ndash;70. <a href=\"https://doi.org/10.1016/j.apenergy.2013.06.057\">https://doi.org/10.1016/j.apenergy.2013.06.057</a>.</li>
<li>Baik, Kyung Don, and Min Soo Kim. &ldquo;Characterization of Nitrogen Gas Crossover through the Membrane in Proton-Exchange Membrane Fuel Cells.&rdquo; <i>International Journal of Hydrogen Energy</i>, 11th International Conference: &ldquo;Hydrogen Materials Science &amp; Chemistry of Carbon Nanomaterials,&rdquo; 36, no. 1 (January 1, 2011): 732&ndash;39. <a href=\"https://doi.org/10.1016/j.ijhydene.2010.09.046\">https://doi.org/10.1016/j.ijhydene.2010.09.046</a>. </li>
<li>Awasthi, A., Keith Scott, and S. Basu. &ldquo;Dynamic Modeling and Simulation of a Proton Exchange Membrane Electrolyzer for Hydrogen Production.&rdquo; <i>International Journal of Hydrogen Energy</i> 36, no. 22 (November 2011): 14779&ndash;86. <a href=\"https://doi.org/10.1016/j.ijhydene.2011.03.045\">https://doi.org/10.1016/j.ijhydene.2011.03.045</a>. </li>
<li>Water Management: EC Kumbur and MM Mench, Fuel Cells &ndash; Proton-Exchange Membrane Fuel Cells, The Pennsylvania State University, University Park, PA, USA</li>
</ol>
</html>", revisions="<html>
Copyright &copy; 2004-2024, MODELON AB <br /> The use of this software component is regulated by the licensing conditions for Modelon Libraries. <br />This copyright notice must, unaltered, accompany all components that are derived from, copied from, <br />or by other means have their origin from any Modelon Library.
</html>"));

end PartialCellTransport_1;
