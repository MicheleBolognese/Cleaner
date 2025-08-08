within PEMFCModel.Reaction.Templates;
partial model BaseReaction "Interface for reaction objects"

  extends PEMFCModel.Reaction.Interfaces.PartialBaseReaction;

  constant .Modelica.Units.SI.Pressure pref = 1e5 "Reference pressure";
  
  Real log10_y_out[N, Medium.nS](start = fill(log10_yout_start, N)) "Logarithm of outlet molar fractions for each control volume";
  Real log10_Q_reac[N, n_reac] "Reaction quotients for each reaction at oulet of each control volume";
  Modelica.Units.SI.MassFraction X_out_real[N, Medium.nS] "Mass-based outlet composition for each control volume after processing";
  Modelica.Units.SI.MassFraction X_out_rescaled[N, Medium.nS] "Mass-based outlet composition for each control volume after rescaling";
  Real check_X_out_real[N](each start = 1) = {sum(X_out_real[j, :]) for j in 1:N} "Check for sum of mass fractions in each control volume";
  Real a[N] = fill(1.0, N)./check_X_out_real "Scaling coefficients for mass fractions if their sum is not 1 in each control volume";
  Real check_X_out_rescaled[N](each start = 1) = {sum(X_out_rescaled[j, :]) for j in 1:N} "Check for sum of mass fractions in each control volume after rescaling";

protected

  Real Sx_reac[n_reac, Medium.nS] "Mass-based stoichiometric matrix";

equation

  // Calculate current outlet composition and reaction quotient
  for i in 1: Medium.nS loop
    
    X_out[:, i] = y_out[:, i]*MMX[i]./MM;
    Sx_reac[:, i] = {S_reac[j, i]*(Medium.MMX[i]/Medium.MMX[comp[j]]) for j in 1:n_reac};

  end for;
  for i in 1:N loop

    log10_Q_reac[i,:] = sum(S_reac[:,j] for j in 1:Medium.nS)*.Modelica.Math.log10(p[i]/pref) + S_reac*log10_y_out[i,:];

  end for; 

algorithm
  
  X_out_real := X_out;
  for j in 1:N loop
        
    for i in 1:Medium.nS loop
   
      if X_out[j, i]-1 > 1e-5 then
     
        assert(false, "Mass fraction of species " + Medium.substanceNames[i] + " in medium " + Medium.mediumName + " in control volume " + String(j) + " is out of the range [0,1]: X = " + String(X_out[j, i]), level = AssertionLevel.warning);
        X_out_real[j, i] := 1 + 1e-5;
     
      end if;
      if X_out[j, i] < -1e-5 then
           
        assert(false, "Mass fraction of species " + Medium.substanceNames[i] + " in medium " + Medium.mediumName + " in control volume " + String(j) + " is out of the range [0,1]: X = " + String(X_out[j, i]), level = AssertionLevel.warning);
        X_out_real[j, i] := -1e-5;

      end if;
  
    end for;
    X_out_rescaled[j, :] := X_out_real[j, :];
    if abs(check_X_out_real[j]-1) > 1e-5 then
   
      assert(false, "Sum of mass fractions in control volume " + String(j) + " is not equal to 1: sum(X_out) = " + String(check_X_out_real[j]), level = AssertionLevel.warning);
      X_out_rescaled[j, :] := a[j]* X_out_real[j, :];

    end if;
    
  end for;

  annotation (Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,
              -100},{100,100}}), graphics={Ellipse(
            extent={{-80,80},{80,-80}},
            lineColor={0,0,0},
            fillPattern=FillPattern.Sphere,
            fillColor={255,128,0})}), Documentation(info="<html>
<h4>BaseReaction</h4>
<p> This is an interface class to general volume reactions. It should only have reactions
and component balances, while it assumes that volume mass and energy balances are handled
outside of the reaction object.
</p>
 
<p> The model where the reaction object is used should propagate parameter values for
initialization etc., and provide the following inputs via modifiers:<br>
<ul>
  <li> T, the reaction temperature, for example wall or in the bulk fluid.</li>
  <li> X_in, mass fractions at the volume inlet. This is used for static conversion models.</li>
  <li> mX_flow, component mass flow balance as the sum of inlet and outlet flows. This is
    used for a dynamic balance in the invariant case.</li>
  <li> reaction parameters, like S_reac, are provided in specific classes inheriting from
    this. In the case of invariant reactions S_reac has to be redeclared as constant.</li>
</ul></p>
 
<p> The inheriting classes provide the actual reaction implementation, and need to give
equations to calculate the following properties:<br>
<ul>
  <li> X_out, mass fractions at the volume outlet.</li>
  <li> Cp_reac, the equivalent heat capacity of the reaction, which is needed for NTU-type HX.</li>
  <li> Q, reaction heat ONLY for the case when heat of formation is not included, default is zero.</li>
  <li> log10_y_out / log10_Q_reac, can not be given here since different variants are
    required for invariant and quasi-dynamic cases. Either one of these should be provided.</li>
</ul></p>
</html>",   revisions="<html>
Copyright &copy; 2004-2024, MODELON AB <br /> The use of this software component is regulated by the licensing conditions for Modelon Libraries. <br />This copyright notice must, unaltered, accompany all components that are derived from, copied from, <br />or by other means have their origin from any Modelon Library.
</html>"));

end BaseReaction;