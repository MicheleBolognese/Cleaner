within PEMFCModel.Reaction.Interfaces;
partial model PartialBaseReaction "Interface for reaction objects"
  
  replaceable package Medium = .FuelCell.Media.PreDefined.IdealGases.NASAReformateLong
  constrainedby .FuelCell.Media.Templates.ReactionGas "Reacting gas properties" annotation(choicesAllMatching);

  parameter Integer N(min = 1) = 1 "Number of control volumes";
  
  input .Modelon.Media.Units.Pressure p[N] "Pressure in each control volume";
  input .Modelon.Media.Units.Temperature T[N] "Reaction temperature in each control volume";
  input .Modelon.Media.Units.MassFraction X_in[N, Medium.nS] "Mass-based inlet composition for each control volume";
  input .Modelon.Media.Units.MassFlowRate mX_flow[N, Medium.nS] "Net mass flow rate of each component through each control volume";
  parameter Real[:, Medium.nS] S_reac = zeros(0, Medium.nS) "Stoichiometric matrix for reactions" annotation(Evaluate=true, Dialog(tab="Advanced"));
  parameter Integer n_reac = size(S_reac, 1) "Number of reactions" annotation(Evaluate=true, Dialog(enable=false, tab="Advanced", group="Predefined"));
  parameter Integer comp[n_reac] = fill(0, n_reac) "Index of reference component for each reaction" annotation(Evaluate=true, Dialog(tab="Advanced"));
  parameter .Modelon.Media.Units.AbsolutePressure pstart = 1e5  "Guessed start value for pressure" annotation(Dialog(group="Initialization"));
  parameter .Modelon.Media.Units.Temperature Tstart = 298.15 "Guessed start value for temperature" annotation(Dialog(group="Initialization"));
  parameter .Modelon.Media.Units.MassFraction Xout_start[Medium.nS] = Medium.reference_X "Guessed start mass-based composition composition at outlet" annotation(Dialog(group="Initialization"));
  parameter .Modelica.Units.SI.Volume V(start = 0.001) = 1 "Total volume";
  parameter Real scale = 1 "Speed factor for quasi-equilibrium reactions" annotation(Dialog(tab="Advanced"));
  parameter String volName = getInstanceName() "Volume-name for better diagnosis" annotation(Dialog(tab="Advanced"));
  
  output .Modelica.Units.SI.HeatFlowRate[N] Q = fill(0, N) "Reaction heat, only if not included in enthalpy";
  output .Modelon.Media.Units.MassFraction[N, Medium.nS] X_out "Mass-based outlet composition for each control volume";

  constant .Modelon.Media.Units.MolarMass MMX[Medium.nS] = Medium.MMX "Molar mass of each medium component" annotation(Dialog(enable = false));
  final parameter Real yout_start[Medium.nS] = Medium.massToMoleFractions(Xout_start, MMX) "Guessed start molar-based composition at outlet" annotation(Dialog(enable = false));
  final parameter Real[Medium.nS] log10_yout_start = {.Modelica.Math.log10(max(.Modelica.Constants.eps, yout_start[i])) for i in 1:Medium.nS} annotation(Dialog(enable = false));
  .Modelon.Media.Units.MoleFraction y_out[N, Medium.nS](start = fill(yout_start, N)) "Outlet molar-based composition for each control volume";
  .Modelon.Media.Units.MolarMass MM[N](start = {1/(Xout_start*(ones(Medium.nS)./MMX)) for i in 1:N}) = {1/(X_out[i,:]*(ones(Medium.nS)./MMX)) for i in 1:N} "Mixture molar mass in each control volume";
  .Modelon.Media.Units.SpecificHeatCapacity[N] Cp_reac "Equivalent specific heat capacity at constant pressure for reactions";

  annotation (Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,
            -100},{100,100}}), graphics={Ellipse(
          extent={{-80,80},{80,-80}},
          lineColor={0,0,0},
          fillPattern=FillPattern.Sphere,
          fillColor={255,128,0})}), Documentation(info="<html>
<h4>PartialBaseReaction</h4>
<p>This is an interface class to general volume reactions (kinetic, equilibrium, invariant). It should only have reactions and component balances, while it assumes that volume mass and energy balances are handled outside of the reaction object. </p>
<p>The model where the reaction object is used should propagate parameter values for initialization etc., and provide the following inputs via modifiers:</p>
<p><ul>
<li>T, the reaction temperature, for example wall or in the bulk fluid. </li>
<li>X_in, mass fractions at the volume inlet. This is used for static conversion models. </li>
<li>dmX_dt_reac, component mass flow balance as the sum of inlet and outlet flows. This is used for a dynamic balance in the invariant case. </li>
<li>reaction parameters, like S_reac, are provided in specific classes inheriting from this. In the case of invariant reactions S_reac has to be redeclared as constant. </li>
</ul></p>
<p>The inheriting classes provide the actual reaction implementation, and need to give equations to calculate the following properties:</p>
<p><ul>
<li>X_out, mass fractions at the volume outlet. </li>
<li>Cp_reac, the equivalent heat capacity of the reaction, which is needed for NTU-type HX. </li>
<li>Q, reaction heat ONLY for the case when heat of formation is not included, default is zero. </li>
<li>log10_y_out / log10_eq_reac, can not be given here since different variants are required for invariant and quasi-dynamic cases. Either one of these should be provided.</li>
</ul></p>
</html>", revisions="<html>
Copyright &copy; 2004-2024, MODELON AB <br /> The use of this software component is regulated by the licensing conditions for Modelon Libraries. <br />This copyright notice must, unaltered, accompany all components that are derived from, copied from, <br />or by other means have their origin from any Modelon Library.
</html>"));

end PartialBaseReaction;
