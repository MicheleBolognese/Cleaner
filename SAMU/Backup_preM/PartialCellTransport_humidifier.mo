within SAMU.Backup_preM;

partial model PartialCellTransport_humidifier "General fuel cell membrane interface with additional transport model"

  extends .SAMU.PartialCell_humidifier;

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

//   parameter Real CF_0_eod = 1.0 "Calibration factor for electro-osmotic drag" annotation(Dialog(tab = "Mass transport", group = "Water transport",enable = true));
//   parameter Real CF_N_eod[N] = fill(CF_0_eod, N) "Calibration factor for electro-osmotic drag (set if non-uniform)" annotation(Dialog(tab = "Mass transport", group = "Water transport",enable = true));    

//   ElectroOsmoticDrag electroOsmoticDrag(
//     redeclare package Medium_an = Medium_an,
//     redeclare package Medium_cath = Medium_cath,
//     final N = N,
//     final n_cell = n_cell,
//     final A_cell = A_cell,
//     final z = z,
//     final rho_dry_m = rho_dry_m,
//     final EW_m = EW_m,
//     final I_cell =  j_ionic*A_cell/N /*fill(I_stack/N, N) */,
//     final p_an_partial = p_an_partial,
//     final p_cath_partial = p_cath_partial,
//     final T_an = T_an,
//     final T_cath = T_cath,
//     final T_cell = T_cell,
//     final iH2O_an = iH2O_an[1],
//     final iH2O_cath = iH2O_cath[1],
//     final lambda = waterContent.lambda,
//     final CF_0 = CF_0_eod,
//     final CF_N = CF_N_eod) annotation (Placement(transformation(extent={{-69.0336325774397,94.04853189230954},{-49.033632577439704,114.04853189230954}},rotation = 0.0,origin = {0.0,0.0})));

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
    final CF_N = CF_N_gasDiffusion) annotation (Placement(transformation(extent={{-68.0,94.0},{-48.0,114.0}},rotation = 0.0,origin = {0.0,0.0})));

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
  
  final .Modelica.Blocks.Interfaces.RealOutput mX_flow_an_transport[N, Medium_an.nS](each unit = "kg/s") "Water/gas transport at dry side (> 0 from dry to wet)";
  final .Modelica.Blocks.Interfaces.RealOutput mX_flow_cath_transport[N, Medium_cath.nS](each unit = "kg/s") "Water/gas transport at wet side (> 0 from wet to dry)";
  final parameter Integer iH2O_an[1] = Medium_an.substanceIndexVector({"H2O"}) "Index of H2O in the dry medium";
  final parameter Integer iH2O_cath[1] = Medium_cath.substanceIndexVector({"H2O"}) "Index of H2O in the wet medium";

protected
  
  final parameter Integer i_required_an[size(diffusiveSpecies, 1)] = Medium_an.substanceIndexVector(diffusiveSpecies) "Indices of required species in the dry medium";
  final parameter Integer i_required_cath[size(diffusiveSpecies, 1)] = Medium_cath.substanceIndexVector(diffusiveSpecies) "Indices of required species in the wet medium";

 initial equation    
  
  for i in 1:size(diffusiveSpecies, 1) loop
  
    assert(i_required_an[i] > 0, "Membrane transport: the selected medium for the dry side does not contain the required substances: see 'required_species'",
    level = AssertionLevel.error);
    assert(i_required_cath[i] > 0, "Membrane transport: the selected medium for the wet side does not contain the required substances: see 'required_species'",
    level = AssertionLevel.error);
  
  end for;

protected
  .Modelica.Units.SI.SpecificEnthalpy h_fg[N] "Enthalpy of vaporization for each node";


equation

  // Mass flow rates of each chemical species to mass flow ports per each stack node, excluding reacting species flows
  for i in 1:N loop
    
    for j in 1:nS_an loop
      
      if j == iH2O_an[1] then

        mX_flow_an_transport[i, j] = waterDiffusion.m_flow[i] /*+ electroOsmoticDrag.m_flow[i] */ /* + waterPermeation.m_flow[i] + gasDiffusion.mX_flow_an[i, j]*/;
      
      else
       
        mX_flow_an_transport[i, j] = gasDiffusion.mX_flow_an[i, j];
      
      end if;
    end for;
    for j in 1:nS_cath loop
      
      if j == iH2O_cath[1] then
        
        mX_flow_cath_transport[i, j] = -waterDiffusion.m_flow[i] /*- electroOsmoticDrag.m_flow[i] */ /* - waterPermeation.m_flow[i] + gasDiffusion.mX_flow_cath[i, j] */;
      
      else
      
        mX_flow_cath_transport[i, j] = gasDiffusion.mX_flow_cath[i, j];
      
      end if;
    
    end for;
  
  end for;
  // HH: Add latent heat transfer for water transport
  for i in 1:N loop
    // Calculate enthalpy of vaporization at membrane temperature
    h_fg[i] = .Modelica.Media.Water.IF97_Utilities.hv_p(satPressure(T_cell[i])) - .Modelica.Media.Water.IF97_Utilities.hl_p(satPressure(T_cell[i]));
    
    // Add latent heat to enthalpy flows (water diffusion > 0 means dry to wet side)
    port_an[i].h_outflow = inStream(port_an[i].h_inflow) - waterDiffusion.m_flow[i] * h_fg / max(1e-6, abs(port_an[i].m_flow));
    port_cath[i].h_outflow = inStream(port_cath[i].h_inflow) + waterDiffusion.m_flow[i] * h_fg / max(1e-6, abs(port_cath[i].m_flow));
  end for;

  annotation (Documentation(info="<html>
<p>Membrane humidifier model with detailed water transport physics. This model has been modified from the original fuel cell transport model to serve as a passive membrane humidifier.</p>
</html>"));

end PartialCellTransport_humidifier;
