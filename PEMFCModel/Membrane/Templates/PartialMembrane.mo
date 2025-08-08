within PEMFCModel.Membrane.Templates;
partial model PartialMembrane "Cell membrane with underlying membrane equations and loss models added, including losses due to contaminants"
  
  extends PEMFCModel.Membrane.Templates.PartialCellTransport;

  parameter .Modelica.Units.SI.Voltage E0_ref = 1.229 "Single cell Nernst's potential at standard conditions" annotation (Dialog(enable = enable_setting,group = "Characteristics and polarization"));
  parameter String an_names[:] = fill("",0) "Required species in the anodic medium"  annotation (Dialog(enable = enable_setting));
  parameter String cath_names[:] = fill("",0) "Required species in the cathodic medium"  annotation (Dialog(enable = enable_setting));
  parameter String contaminantsSpecies[:] = fill("", 0) "Species affecting voltage loss due to contaminants" annotation(Dialog(enable = enable_setting,group = "Characteristics and polarization"));
  parameter Real n_e = 2.0 "Number of exchanged electrons per unit amount of reference chemical species" annotation (Evaluate = true, Dialog(enable = enable_setting,tab = "Reactions"));
  parameter Real S_reac_an[nS_an] = Medium_an.stoichiometry(an_names, fill(1, size(an_names,1))) "Stoichiometry vector for reactions at the anode side: - if consumed, + if generated" annotation (Evaluate = true, Dialog(tab = "Reactions", enable = enable_setting));
  parameter Real S_reac_cath[nS_cath] = Medium_cath.stoichiometry(cath_names, fill(1, size(cath_names,1))) "Stoichiometry vector for reactions at the cathode side: - if consumed, + if generated" annotation(Evaluate = true, Dialog(tab = "Reactions",enable = enable_setting));
   
  .Modelica.Units.SI.Voltage E0_cell[N] "Single cell Nernst's potential in each node at operating conditions";
  .Modelica.Units.SI.Voltage E_act_cell[N] "Single cell activation loss in each node";
  .Modelica.Units.SI.Voltage E_conc_cell[N] "Single cell concentration loss in each node";
  .Modelica.Units.SI.Voltage E_ohm_cell[N] "Single cell ohmic loss in each node";
  .Modelica.Units.SI.Voltage E_cont_cell[N] "Single cell voltage loss due to contaminants in each node";
 
  Real eff_cell "Cell efficiency in terms of output power/energy content of consumed hydrogen";
  Real eff_volt_cell "Cell efficiency in terms of output voltage/ideal voltage";

  .Modelica.Units.SI.MassFlowRate mX_flow_an_reac[N, nS_an] "Mass flow rates of each reacting species only from anode channel in each node for the stack (> 0 from anode channel to membrane)";
  .Modelica.Units.SI.MassFlowRate mX_flow_cath_reac[N, nS_cath] "Mass flow rates of each reacting species only from cathode channel in each node for the stack (> 0 from cathode channel to membrane)";

  .Modelica.Electrical.Analog.Sources.SignalVoltage sourceV[N] annotation (Placement(transformation(
        origin={0.0,0.0},
        extent={{-0.11146927037651366,-40.14719423510464},{19.888530729623486,-60.14719423510464}},
        rotation=0.0)));
  .Modelica.Blocks.Sources.RealExpression Voc_stack[N](y = E0_cell*n_cell) annotation (Placement(transformation(
        extent={{34.62471514085864,-87.68058263722398},{18.624715140858644,-75.68058263722398}},
        rotation=0.0,
        origin={0.0,0.0})));

  replaceable model ActivationLoss = .PEMFCModel.Membrane.Losses.ActivationLoss.ZeroLoss constrainedby
    .PEMFCModel.Membrane.Losses.Templates.BaseActivationLoss(final enableInternal = false) annotation (choicesAllMatching, Dialog(tab = "Voltage losses",enable = true));

  replaceable model OhmicLoss = .PEMFCModel.Membrane.Losses.OhmicLoss.ZeroLoss constrainedby         
    .PEMFCModel.Membrane.Losses.Templates.BaseOhmicLoss(final enableInternal = false) annotation (choicesAllMatching, Dialog(tab = "Voltage losses",enable = true));

  replaceable model ConcentrationLoss = .PEMFCModel.Membrane.Losses.ConcentrationLoss.ZeroLoss constrainedby
    .PEMFCModel.Membrane.Losses.Templates.BaseConcentrationLoss(final enableInternal = false)  annotation (choicesAllMatching, Dialog(tab = "Voltage losses",enable = true));
  
  replaceable model ContaminantsLoss = .PEMFCModel.Membrane.Losses.ContaminantsLoss.ZeroLoss constrainedby .PEMFCModel.Membrane.Losses.Templates.BaseContaminantsLoss(final enableInternal = false) annotation (choicesAllMatching, Dialog(tab = "Voltage losses",enable = true));
  
   ActivationLoss activationLoss(
    final N = N,
    final n_cell = n_cell,
    final j_ionic = j_ionic) "Activation loss" annotation (Placement(transformation(
        extent={{-10.0,-10.0},{10.0,10.0}},
        rotation=-90.0,
        origin={41.432909514482134,33.84511874268364})));

  OhmicLoss ohmicLoss(
    final N = N,
    final n_cell = n_cell,
    final A_cell = A_cell) "Ohmic loss" annotation (Placement(transformation(
        extent={{-10.0,10.0},{10.0,-10.0}},
        rotation=-90.0,
        origin={41.387649370066804,3.161408846947891})));

 ConcentrationLoss concentrationLoss(
    final N = N,
    final n_cell = n_cell,
    final j_ionic = j_ionic) "Concentration loss" annotation (Placement(transformation(
        extent={{-10.0,-10.000000000000002},{10.0,10.000000000000002}},
        rotation=-90.0,
        origin={40.94142259414224,-28.60669456066948})));

  ContaminantsLoss contaminantsLoss(
    redeclare package Medium_an = Medium_an,
    contaminantsSpecies = contaminantsSpecies,
    final N = N,
    final n_cell = n_cell,
    final T_cell = T_cell,
    final p_an_partial = p_an_partial,
    //final p_cath_partial = p_cath_partial,
    final y_an = y_an) "Voltage loss due to contaminants" annotation (Placement(transformation(extent = {{-0.04943406719663912,40.180403110838114},{19.95056593280336,60.180403110838114}},origin = {0.0,0.0},rotation = 0.0)));
 
  .Modelica.Units.SI.MassFlowRate checkMassBalance = abs(sum(port_an.m_flow + port_cath.m_flow));
  .Modelica.Units.SI.Power dUdt = sum(C_cell*n_cell/N*dTdt) "Time-derivative of membrane internal energy";  
  .Modelica.Units.SI.Power checkEnergyBalance = abs(dUdt - sum(port_an.H_flow + port_cath.H_flow) + P_stack - (Q_wall_an_stack + Q_wall_cath_stack + Q_wall_stack));
    
protected
  
  parameter Integer iAn[size(an_names, 1)] = Medium_an.substanceIndexVector(an_names) "Indices of required species in the anodic medium";
  parameter Integer iCath[size(cath_names, 1)] = Medium_cath.substanceIndexVector(cath_names) "Indices of required cathode species in the cathodic medium";

initial equation

  for i in 1:size(an_names,1) loop 
    
    assert(iAn[i] > 0, "Membrane: the selected medium for the anode side does not contain the required substances: see 'an_names'", level = AssertionLevel.error);

  end for;
  for i in 1:size(cath_names,1) loop
    
    assert(iCath[i] > 0, "Membrane: the selected media for the anode side does not contain the required substances: see 'cath_names'", level = AssertionLevel.error);
  
  end for;

equation
 
  //Single cell voltage losses
  E_ohm_cell = ohmicLoss.E_loss_cell;
  E_act_cell = activationLoss.E_loss_cell;
  E_conc_cell = concentrationLoss.E_loss_cell;
  E_cont_cell = contaminantsLoss.E_loss_cell;
  
  //Mass flow rates of reacting chemical species to mass transfer ports
  for i in 1:N loop
    
    mX_flow_an_reac[i, :] = j_ionic[i]/(n_e*.FuelCell.Internal.Units.F)*A_cell/N*n_cell*(-S_reac_an .* Medium_an.MMX);
    mX_flow_cath_reac[i, :] = j_ionic[i]/(n_e*.FuelCell.Internal.Units.F)*A_cell/N*n_cell*(-S_reac_cath .* Medium_cath.MMX);
    for j in 1:nS_an loop
      
      port_an[i].mX_flow[j] = mX_flow_an_reac[i, j] + mX_flow_an_transport[i, j];
    
    end for;
    port_an[i].m_flow = sum(mX_flow_an_reac[i, :]) + sum(mX_flow_an_transport[i, :]);
    port_an[i].H_flow = sum(anode[i].h_component .* mX_flow_an_reac[i, :]) + sum(anode[i].h_component .* mX_flow_an_transport[i,:]);
    for j in 1:nS_cath loop
      
      port_cath[i].mX_flow[j] = mX_flow_cath_reac[i, j] + mX_flow_cath_transport[i, j];
    
    end for;
    port_cath[i].m_flow = sum(mX_flow_cath_reac[i, :]) + sum(mX_flow_cath_transport[i, :]);
    port_cath[i].H_flow = sum(cathode[i].h_component .*mX_flow_cath_reac[i, :]) + sum(cathode[i].h_component .*mX_flow_cath_transport[i, :]);
  
  end for;
   
  ohmicLoss.I = j_ionic*A_cell/N;
 
  n_cell*Q_cell = port_an.H_flow + port_cath.H_flow - n_cell*P_cell;
  P_cell = V_cell*ohmicLoss.I;
  Q_stack = sum(n_cell*Q_cell);  

  eff_cell = max(0, min(1, P_stack/max(.Modelica.Constants.eps, abs(sum(port_an.H_flow + port_cath.H_flow)))));
  eff_volt_cell = max(0, min(1, V_cell/E0_ref));

  for i in 1:N loop
    
    connect(pin_n, contaminantsLoss.pin_p[i]) annotation (Line(points={{-90,50},{0,50}}, color={0,0,255}));
    connect(pin_p, sourceV[i].p) annotation (Line(points={{0,-50},{-90,-50}}, color={0,0,255}));
  
  end for;
  connect(Voc_stack.y,sourceV.v) annotation(Line(points = {{17.82471514085864,-81.68058263722398},{9.888530729623486,-81.68058263722398},{9.888530729623486,-62.14719423510464}},color = {0,0,127}));
  connect(contaminantsLoss.pin_n,activationLoss.pin_p) annotation (Line(points = {{19.95056593280336,50.180403110838114},{41.432909514482134,50.180403110838114},{41.432909514482134,43.84511874268364}},color = {0,0,255}));
  connect(activationLoss.pin_n, ohmicLoss.pin_p) annotation (Line(points = {{41.432909514482134,23.84511874268364},{41.387649370066804,23.84511874268364},{41.387649370066804,13.161408846947891}},color = {0,0,255}));
  connect(ohmicLoss.pin_n, concentrationLoss.pin_p) annotation (Line(points={{41.387649370066804,-6.838591153052109},{40.94142259414224,-6.838591153052109},{40.94142259414224,-18.606694560669478}}, color={0,0,255}));
  connect(concentrationLoss.pin_n, sourceV.n) annotation (Line(points={{40.94142259414224,-38.60669456066948},{40.94142259414224,-50.14719423510464},{19.888530729623486,-50.14719423510464}}, color={0,0,255}));
   
  annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)),
    Documentation(info="<html>
<h4>PartialMembrane</h4>
<p>Partial cell membrane model with mass and energy balances defined based on variable selections.</p>
<p>To create a complete model:</p>
<ul>
<li>Proper replaceable models should be selected and all parameters properly tuned (e.g., for PEMFC vs SOFC).</li>
<li>One additional equation is required to define the open circuit voltage (e.g, nernst equation - E0) in the extending models</li>
</ul>
</html>", revisions="<html>
Copyright &copy; 2004-2024, MODELON AB <br /> The use of this software component is regulated by the licensing conditions for Modelon Libraries. <br />This copyright notice must, unaltered, accompany all components that are derived from, copied from, <br />or by other means have their origin from any Modelon Library.
</html>"));

end PartialMembrane;
