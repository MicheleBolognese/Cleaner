within PEMFCModel.Stacks.Interfaces;
partial model PartialSubStack "Generic substack interface"

  replaceable package Medium_an = .FuelCell.Media.PreDefined.IdealGases.NASAReformateLong constrainedby 
    .FuelCell.Media.Templates.ReactionGas "Anodic medium" annotation (choicesAllMatching, Dialog(enable_setting, tab = "Substack"));

  replaceable package Medium_cath = .FuelCell.Media.PreDefined.IdealGases.NASAMoistAir constrainedby 
    .FuelCell.Media.Templates.ReactionGas "Cathodic medium" annotation (choicesAllMatching,Dialog(enable_setting, tab = "Substack"));

  replaceable model Membrane = PEMFCModel.Membrane.Templates.PartialMembrane(final enable_setting = false) constrainedby PEMFCModel.Membrane.Templates.PartialCell(
    final enable_setting = false) "Cell membrane model to be used" annotation (choicesAllMatching, Dialog(tab = "Cell model",group = "Membrane characteristics"));

  parameter Boolean enable_setting = true "If true, selectability of parameters enabled" annotation(Dialog( tab = "Substack"));

  parameter Integer N(min = 1) = 2 "Number of along-the-channel discretization nodes" annotation(Dialog(enable = enable_setting,tab = "Substack"));
  constant Integer nS_an = Medium_an.nS "Number of chemical species in anodic medium";
  constant Integer nS_cath = Medium_cath.nS "Number of chemical species in cathodic medium";
  parameter Integer n_cell = 50 "Number of series-connected cells in substack" annotation (Dialog(enable = enable_setting,tab = "Substack"));
  parameter .Modelica.Units.SI.Area A_cell = 180e-4 "Cell active area" annotation (Dialog(enable = enable_setting,tab = "Substack"));
  parameter .Modelica.Units.SI.Mass M_stack = 30 "Substack mass" annotation (Dialog(enable = enable_setting,tab = "Substack"));
  parameter .Modelica.Units.SI.SpecificHeatCapacity c_stack = 500 "Specific heat capacity of substack material" annotation (Dialog(enable = enable_setting,tab = "Substack")); 

  parameter Real z = 20e-6 "Membrane thickness" annotation (Dialog(tab = "Cell model", group = "Membrane characteristics", enable = enable_setting));
  parameter Real rho_dry_m = 1980 "Dry membrane density" annotation (Dialog(tab = "Cell model", group = "Membrane characteristics", enable = enable_setting));
  parameter Real EW_m = 1.1 "Membrane equivalent weight" annotation (Dialog(tab = "Cell model", group = "Membrane characteristics", enable = enable_setting));
  parameter Real E0_ref = 1.229 "Single cell Nernst's potential at standard conditions" annotation (Dialog(tab = "Cell model", group = "Characteristics and polarization", enable = enable_setting));
  parameter String contaminantsSpecies[:] = fill("",0) "Species affecting voltage loss due to contaminants" annotation (Dialog(tab = "Cell model", group = "Characteristics and polarization", enable = enable_setting));
  parameter String diffusiveSpecies[:] = fill("",0) "Species diffusing through membrane" annotation (Dialog(tab = "Cell model", group = "Membrane characteristics", enable = enable_setting));
  
  parameter Boolean includeCellConduction = false "If true, along-the-channel thermal conduction included" annotation(Dialog(tab = "Cell model",group = "Along-the-channel thermal conduction",enable = enable_setting));
  parameter .Modelica.Units.SI.ThermalConductivity lambda_cell = 20 "Cell internal thermal conductivity" annotation (Dialog(enable = includeCellConduction, group = "Along-the-channel thermal conduction", tab = "Cell model"));
  parameter .Modelica.Units.SI.Area A_crosssection_cell = 18e-3*2e-3 "Cross section area of a single cell" annotation (Dialog(enable = includeCellConduction, group = "Along-the-channel thermal conduction", tab = "Cell model"));
  parameter .Modelica.Units.SI.Length length_cell = 0.15 "Cell length" annotation (Dialog(enable = includeCellConduction, tab = "Cell model",group = "Along-the-channel thermal conduction"));
  final parameter .Modelica.Units.SI.ThermalConductance G_cell = lambda_cell*A_crosssection_cell/(length_cell/N) "Cell thermal conductance" annotation(Dialog(tab = "Cell model", group = "Along-the-channel thermal conduction", enable = false));

  parameter Boolean addProxToAnode = false "Add prox reactor before inlet to anode channel" annotation (Dialog(tab = "Advanced", group = "Reactions", enable = enable_setting));

  .Modelica.Units.SI.Temperature T_stack[N] "Substack temperature in each node";
  .Modelica.Units.SI.Voltage V_stack "Substack voltage";
  .Modelica.Units.SI.Power P_stack "Substack power";
  
  // Stoichiometry & utilization
  input Modelica.Units.SI.MassFraction X_feed_an[nS_an] "Mass-based composition in feed connector at anodic side";
  .Modelica.Units.SI.MoleFraction anode_stoich "Anode stoichiometry";
  input Modelica.Units.SI.MassFraction X_feed_cath[nS_cath] "Mass-based composition in feed connector at cathodic side";
  .Modelica.Units.SI.MoleFraction cathode_stoich "Cathode stoichiometry";

  Membrane cell(
    N = N,
    A_cell = A_cell,
    n_cell = n_cell,
    C_cell = (M_stack*c_stack)/n_cell,
    redeclare package Medium_an = Medium_an,
    redeclare package Medium_cath = Medium_cath,
    z = z,
    rho_dry_m = rho_dry_m,
    EW_m = EW_m,
    E0_ref = E0_ref,
    contaminantsSpecies = contaminantsSpecies,
    diffusiveSpecies = diffusiveSpecies,T_from_h = false) annotation (Placement(transformation(extent={{-20.310787158962256,-19.49372383321096},{19.689212841037744,20.50627616678904}},
                       rotation=0.0,origin = {0.0,0.0})));

  .Modelica.Electrical.Analog.Interfaces.PositivePin pin_p annotation (Placement(transformation(extent={{-10.0,-78.07112957488314},{10.0,-58.071129574883145}}, rotation=
            0.0,origin = {0.0,0.0}), iconTransformation(extent={{-10,-76},{10,-56}})));

  .Modelica.Electrical.Analog.Interfaces.NegativePin pin_n annotation (Placement(transformation(extent={{-10.0,56.0},{10.0,76.0}}, rotation=0.0,origin = {0.0,0.0}),
        iconTransformation(extent={{-10,56},{10,76}}))) ;

  .FuelCell.Interfaces.GasVolumePort feed_anode(redeclare replaceable package Medium = Medium_an) annotation (Placement(transformation(extent={{-100,30},{-80,50}}, rotation=
            0), iconTransformation(extent={{-100,30},{-80,50}})));

  .FuelCell.Interfaces.GasFlowPort drain_anode(redeclare replaceable package Medium = Medium_an) annotation (Placement(transformation(extent={{80,30},{100,50}}, rotation=0),
        iconTransformation(extent={{80,30},{100,50}})));

  .FuelCell.Interfaces.GasVolumePort feed_cathode(redeclare replaceable package Medium = Medium_cath) annotation (Placement(transformation(extent={{-98,-50},{-78,-30}},
          rotation=0), iconTransformation(extent={{-100,-50},{-80,-30}})));

  .FuelCell.Interfaces.GasFlowPort drain_cathode(redeclare replaceable package Medium = Medium_cath) annotation (Placement(transformation(extent={{80,-50},{100,-30}}, rotation=
            0), iconTransformation(extent={{80,-50},{100,-30}})));

  .Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a wall[N](each T(nominal = 500)) annotation (Placement(transformation(extent={{30,-10},
            {50,10}},          rotation=0), iconTransformation(extent={{30,-10},
            {50,10}})));

  /* Thermal conductors along flow direction */
  .Modelica.Thermal.HeatTransfer.Components.ThermalConductor[N-1] conduction(each G = G_cell*n_cell) if includeCellConduction and N > 1 "Thermal conduction along flow direction";

equation

  for i in 1:N-1 loop

    connect(cell.wall[i], conduction[i].port_a);
    connect(conduction[i].port_b, cell.wall[i+1]);

  end for;

  annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
            -100},{100,100}})),
                       Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,
            -100},{100,100}}), graphics={
        Rectangle(
          extent={{-68,24},{70,-24}},
          lineColor={0,0,0},
          fillColor={135,135,135},
          fillPattern=FillPattern.CrossDiag),
        Rectangle(
          extent={{-80,56},{80,24}},
          lineColor={0,0,0},
          fillPattern=FillPattern.HorizontalCylinder,
          fillColor={215,215,215}),
        Rectangle(
          extent={{-80,-24},{80,-56}},
          lineColor={0,0,0},
          fillPattern=FillPattern.HorizontalCylinder,
          fillColor={215,215,215}),
        Text(
          extent={{-98,20},{100,8}},
          lineColor={255,255,255},
          textString="%name")}),
    Documentation(info="<html>
<h4>PartialSubStack</h4>
<p>Interface model to be used as a basis for the Substack template. The model is an interface class which means that the equations needs to be finalized before the model can be used. A finalized templates that inherits from the class is found in FuelCell.Stacks.Templates.
The model contains the following components:</p>
<p><ul>
<li>Feed and drain connectors of reaction gas type (idealgas)</li>
<li>Electrical pins (positive and negative)</li>
<li>Replaceable partial cell membrane model (used as placeholder) - needs to be replaced before simulation</li>
<li>Heat port connector (wall), discretized in the flow direction to describe the temperature profile in the flow direction of the membrane </li>
</ul></p>
<p>The model also contains parameters and variables common for all Substack templates (stack voltage, stack temperature etc), parameters for model discretization and anode and cathode media models.</p>
<p></p>
</html>", revisions="<html>
Copyright &copy; 2004-2024, MODELON AB <br /> The use of this software component is regulated by the licensing conditions for Modelon Libraries. <br />This copyright notice must, unaltered, accompany all components that are derived from, copied from, <br />or by other means have their origin from any Modelon Library.
</html>"));

end PartialSubStack;
