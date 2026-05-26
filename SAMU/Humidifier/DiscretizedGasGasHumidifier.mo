within SAMU.Humidifier;

model DiscretizedGasGasHumidifier
  "Model for a discretized gas-gas humidifier"
  extends .FuelCell.HeatExchangers.Interfaces.DistributedGasGasHumidifier(
    wallThickness=t_mem,
    redeclare replaceable model HeatTransfer_sec =
        .Modelon.ThermoFluid.FlowChannels.HeatTransfer.SinglePhase.DittusBoelterAdjustable,
    redeclare replaceable model Friction_sec =
        .Modelon.ThermoFluid.FlowChannels.PipeResistances.SinglePhase.QuadraticOperatingPointLoss,
    redeclare replaceable model HeatTransfer_prim =
        .Modelon.ThermoFluid.FlowChannels.HeatTransfer.SinglePhase.DittusBoelterAdjustable,
    redeclare replaceable model Friction_prim =
        .Modelon.ThermoFluid.FlowChannels.PipeResistances.SinglePhase.QuadraticOperatingPointLoss,
    redeclare
      .Modelon.ThermoFluid.HeatExchangers.Records.Summary.Base_pinchlmtd.BaseSummary
      summary(
      m_flow=m_flow_prim,
      m_flow_sec=m_flow_sec,
      Q_flow=sum(shell.wall.Q_flow),
      T_in=noEvent(if m_flow_prim >= 0 then shell.summary.T_in else shell.summary.T_out),
      T_out=noEvent(if m_flow_prim >= 0 then shell.summary.T_out else shell.summary.T_in),
      h_in=noEvent(if m_flow_prim >= 0 then shell.summary.h_in else shell.summary.h_out),
      h_out=noEvent(if m_flow_prim >= 0 then shell.summary.h_out else shell.summary.h_in),
      dp=shell.summary.dp,
      T_sec_in=noEvent(if m_flow_sec >= 0 then tube.summary.T_in else tube.summary.T_out),
      T_sec_out=noEvent(if m_flow_sec >= 0 then tube.summary.T_out else tube.summary.T_in),
      h_sec_in=noEvent(if m_flow_sec >= 0 then tube.summary.h_in else tube.summary.h_out),
      h_sec_out=noEvent(if m_flow_sec >= 0 then tube.summary.h_out else tube.summary.h_in),
      dp_sec=tube.summary.dp,
      p_in=portA_prim.p,
      p_out=portB_prim.p,
      p_sec_in=portA_sec.p,
      p_sec_out=portB_sec.p,
      Q_flow_sec=sum(tube.wall.Q_flow)));

  /* Local variables */
  .Modelica.Units.SI.MassFlowRate m_flow_prim
    "Mass flowrate primary side, positive flow from portA_prim to portB_prim";
  .Modelica.Units.SI.MassFlowRate m_flow_sec
    "Mass flowrate primary side, positive flow from portA_sec to portB_sec";
  .Modelica.Units.SI.HeatFlowRate Q_tot "Total heat flow rate out from primary side";

  .FuelCell.Pipes.FlowChannel shell(
    redeclare package Medium = PrimaryMedium,
    A=fill(.Modelica.Constants.pi*D_h*D_h/4, n),
    A_heat=fill(A_mem, n),
    useHeatTransfer=useHeatTransfer,
    initOpt=initOpt,
    p_start_in=init_prim.p_in,
    p_start_out=init_prim.p_out,
    p_start=p_start,
    initFromEnthalpy=true,
    h_start_in=init_prim.h_in_used,
    h_start_out=init_prim.h_out_used,
    h_start=h_start,
    X_start=init_prim.X,
    m_flow_start=init_prim.m_flow,
    L=fill(L_mem/n, n),
    Dhyd=fill(D_h, n),
    positiveFlow=positiveFlow_prim,
    from_dp=from_dp_prim,
    dp_smooth=dp_smooth,
    mflow_smooth=mflow_smooth,
    generateEventForReversal=generateEventForReversal_prim,
    C_start=init_prim.C,
    redeclare model Friction = Friction_prim,
    redeclare model HeatTransfer = HeatTransfer_prim,
    CF_PressureLoss=CF_PrimarySidePressureLoss,
    CF_HeatTransfer=CF_PrimarySideHeatTransfer,
    //H_flow = h_water_prim .* m_trans_shell,
    mX_flow=transpose({if i == H2O_prim then -m_trans_shell else zeros(n) for i in 1:
        nX_shell})) "Humidifier shell side" annotation (Dialog(tab="Sub-components",
        group="Flow channels"), Placement(transformation(
        extent={{-20,-20},{20,20}},
        rotation=180,
        origin={-8,60})));

  .FuelCell.Pipes.FlowChannel tube(
    redeclare package Medium = SecondaryMedium,
    n_channels=fill(n_tubes, n),
    A=fill(.Modelica.Constants.pi*D_mem*D_mem/4, n),
    A_heat=fill(A_mem, n),
    useHeatTransfer=useHeatTransfer,
    initOpt=initOpt_sec,
    p_start_out=init_sec.p_out_used,
    p_start=p_sec_start,
    initFromEnthalpy=true,
    h_start_in=init_sec.h_in_used,
    h_start_out=init_sec.h_out_used,
    h_start=h_sec_start,
    X_start=init_sec.X,
    m_flow_start=init_sec.m_flow,
    p_start_in=init_sec.p_in_used,
    L=fill(L_mem/n, n),
    Dhyd=fill(D_mem, n),
    positiveFlow=positiveFlow_sec,
    from_dp=from_dp_sec,
    dp_smooth=dp_smooth,
    mflow_smooth=mflow_smooth,
    generateEventForReversal=generateEventForReversal_sec,
    C_start=init_sec.C,
    redeclare model Friction = Friction_sec,
    redeclare model HeatTransfer = HeatTransfer_sec,
    CF_PressureLoss=CF_SecondarySidePressureLoss,
    CF_HeatTransfer=CF_SecondarySideHeatTransfer,
    //H_flow = h_water_sec .* m_trans,
    mX_flow=transpose({if i == H2O_sec then m_trans else zeros(n) for i in 1:nX_tube}))
    "Humidifier tube side" annotation (Dialog(tab="Sub-components", group=
          "Flow channels"), Placement(transformation(extent={{-28,-80},{12,-40}})));

  .Modelon.ThermoFluid.Solids.DynamicWall wall(
    props(each Cp=wallMaterial.c, Rw=wallThickness./wallMaterial.lambda./
          wallCrossSecArea),
    n=n,
    steadyStateInit=WallSteadyStateInit,
    n_channels=n_channels_wall,
    m=m,
    T0=T0_wall,
    massLessWall=massLessWall,
    includeThermalResistance=includeThermalResistance,
    surfaceState = false) "Wall model definition"
    annotation (Dialog(tab="Sub-components", group="Wall"), Placement(
        transformation(extent={{-30.0,-20.0},{10.0,20.0}},rotation = 0.0,origin = {0.0,0.0})));

  replaceable record WallMaterial =
      .Modelon.Thermal.MaterialProperties.PropertyData.ConstantProperties.SolidMaterial
      (
      c=4180,  
      rho= 1968,
      lambda=0.12) constrainedby
    .Modelon.Thermal.MaterialProperties.PropertyData.ConstantProperties.SolidMaterial
    "Geometry and material parameters"     annotation(choicesAllMatching, Dialog(tab="Wall",group="Wall geometry"));

  WallMaterial wallMaterial;

  parameter Boolean useHeatTransfer=false
    "Consider heat transfer effects (not used when connected to membrane)"
    annotation (Dialog(tab="Correlations", group="Heat transfer"));

  /* Diffusion correlation */

  constant Integer H2O_prim=PrimaryMedium.substanceIndex("H2O");
  constant Integer H2O_sec=SecondaryMedium.substanceIndex("H2O");
  constant Integer nX_shell = PrimaryMedium.nS
    "total number of mass fractions for shell side";
  constant Integer nX_tube = SecondaryMedium.nS
    "total number of mass fractions for tube side";

  package WaterMedium = .Modelon.Media.PreDefined.TwoPhase.WaterIF97
    "water medium model";

  .Modelica.Units.SI.MassFlowRate m_water=sum(m_trans) "Water diffusion mass flow rate";
  .Modelica.Units.SI.MassFlowRate[n] m_trans "Mass flow rate of steam across membrane per segment";

  .Modelon.Media.Units.MoleFraction[n,nX_shell] y_prim
     "primary inlet concentration";
  .Modelon.Media.Units.MoleFraction[n,nX_tube] y_sec
     "secondary exit concentration";

  .Modelica.Units.SI.DimensionlessRatio a_shell[n] "Water activity of shell side";
  .Modelica.Units.SI.DimensionlessRatio a_tube[n] "Water activity of tube side";
  .Modelica.Units.SI.DimensionlessRatio a_mem[n] "Water activity of membrane side";

  Real lambda_mem[n] "Mean water content in membrane";
  Real lambda_shell[n] "Water content in shell side";
  Real lambda_tube[n] "Water content in tube side";

  .Modelica.Units.SI.Temperature T_mem[n] "Membrane temperature";

  PrimaryMedium.ThermodynamicState state_prim[n]
     "thermodynamic state on the primary side";
  SecondaryMedium.ThermodynamicState state_sec[n]
     "thermodynamic state on the secondary side";

  .Modelica.Units.SI.DiffusionCoefficient Dw[n] "Membrane diffusion coefficient";

//protected
  .Modelica.Units.SI.MassFlowRate[n] m_trans_shell "Mass transfer rate for shell";
  Real D[n] "Empirical constant";

  .Modelica.Units.SI.Pressure pV_tube[n] "Partial pressure of water vapor tube side";
  .Modelica.Units.SI.Pressure pV_shell[n] "Partial pressure of water vapor shell side";
  .Modelica.Units.SI.Pressure psat_tube[n] "Water saturation pressure tube side";
  .Modelica.Units.SI.Pressure psat_shell[n] "Water saturation pressure shell side";

  .Modelica.Units.SI.MolarDensity C_shell[n] "Water mass concentration shell side";
  .Modelica.Units.SI.MolarDensity C_tube[n] "Water mass concentration tube side";
  .Modelica.Units.SI.MassConcentration C_w[n] "Overall water mass concentration";

  .Modelica.Units.SI.MolarMass Mv=18.01528e-3 "Molar mass of vapor";

  .Modelon.Media.Units.EnthalpyFlowRate[n] h_water_prim
     "enthalpy of gaseous water";
  .Modelon.Media.Units.EnthalpyFlowRate[n] h_water_sec
     "enthalpy of gaseous water";
//   Modelon.Media.Units.EnthalpyFlowRate[n] h_water_prim_g
//      "enthalpy of gaseous water";
//   Modelon.Media.Units.EnthalpyFlowRate[n] h_water_sec_g
//      "enthalpy of gaseous water";
//   Modelon.Media.Units.EnthalpyFlowRate[n] h_water_prim_l
//      "enthalpy of gaseous water";
//   Modelon.Media.Units.EnthalpyFlowRate[n] h_water_sec_l
//      "enthalpy of gaseous water";
initial equation
  assert(
    H2O_prim > 0,
    "Humidifier: PrimaryMedium does not contain the required substance: H2O",
    level=AssertionLevel.error);
  assert(
    H2O_sec > 0,
    "Humidifier: SecondaryMedium does not contain the required substance: H2O",
    level=AssertionLevel.error);

equation
  m_flow_prim =shell.summary.m_flow;
  m_flow_sec =tube.summary.m_flow;
  Q_tot =-1*sum(shell.wall.Q_flow);

  state_sec = tube.channel.volume.state;

 for i in 1:n loop

//     h_water_prim[i] = PrimaryMedium.specificEnthalpy_index(state_prim[i],H2O_prim)*(state_prim[i].X[H2O_prim]);
//    h_water_sec[i] = SecondaryMedium.specificEnthalpy_index(state_sec[i],H2O_sec)*(state_sec[i].X[H2O_sec]);

    h_water_prim[i] = (PrimaryMedium.enthalpyOfCondensingGas(state_prim[i].T) + PrimaryMedium.enthalpyOfLiquid(state_prim[i].T))*(state_prim[i].X[H2O_prim]);
    h_water_sec[i] = (SecondaryMedium.enthalpyOfCondensingGas(state_sec[i].T) + SecondaryMedium.enthalpyOfLiquid(state_sec[i].T))*(state_prim[i].X[H2O_prim]);   // Mod for consistency

     // h_water_prim_g[i] = PrimaryMedium.enthalpyOfCondensingGas(state_prim[i].T);
     // h_water_sec_g[i] = SecondaryMedium.enthalpyOfCondensingGas(state_sec[i].T); 
//
//     h_water_prim_l[i] = PrimaryMedium.enthalpyOfLiquid(state_prim[i].T);
//     h_water_sec_l[i] = SecondaryMedium.enthalpyOfLiquid(state_sec[i].T);

    y_prim[i,:] = PrimaryMedium.massToMoleFractions(state_prim[i].X,PrimaryMedium.MMX);
    y_sec[i,:] = SecondaryMedium.massToMoleFractions(state_sec[i].X,SecondaryMedium.MMX);

    pV_shell[i] = y_prim[i, H2O_prim]*state_prim[i].p;
    pV_tube[i] = y_sec[i, H2O_sec]*state_sec[i].p;

    psat_shell[i] = min(WaterMedium.saturationPressure_TX(max(shell.channel.volume[i].T, 273.16)), 0.999*state_prim[i].p);
    psat_tube[i] = min(WaterMedium.saturationPressure_TX(max(tube.channel.volume[i].T, 273.16)), 0.999*state_sec[i].p);

    a_shell[i] = pV_shell[i]/psat_shell[i];
    a_tube[i] = pV_tube[i]/psat_tube[i];
    a_mem[i] = (a_shell[i] + a_tube[i])/2;

    if a_shell[i] >= 0 and a_shell[i] <= 1 then
      lambda_shell[i] =0.043 + 17.81*a_shell[i] - 39.85*a_shell[i]^2 + 36.0*a_shell[i]^3;
    elseif a_shell[i] > 1 and a_shell[i] <= 3 then
      lambda_shell[i] =14 + 1.4*(a_shell[i] - 1);
    else
      lambda_shell[i] = 16.8;
    end if;

    if a_tube[i] >= 0 and a_tube[i] <= 1 then
      lambda_tube[i] =0.043 + 17.81*a_tube[i] - 39.85*a_tube[i]^2 + 36.0*a_tube[i]^3;
    elseif a_tube[i] > 1 and a_tube[i] <= 3 then
      lambda_tube[i] =14 + 1.4*(a_tube[i] - 1);
    else
      lambda_tube[i] = 16.8;
    end if;

    if a_mem[i] >= 0 and a_mem[i] <= 1 then
      lambda_mem[i] =0.043 + 17.81*a_mem[i] - 39.85*a_mem[i]^2 + 36.0*a_mem[i]^3;
    elseif a_mem[i] > 1 and a_mem[i] <= 3 then
      lambda_mem[i] =14 + 1.4*(a_mem[i] - 1);
    else
      lambda_mem[i] = 16.8;
    end if;

    if lambda_mem[i] < 2 then
      D[i] = 1e-6;
    elseif lambda_mem[i]>=2 and lambda_mem[i]<=3 then
      D[i] = 1e-6*(1 + 2*(lambda_mem[i] - 2));
    elseif lambda_mem[i]> 3 and lambda_mem[i]<4.5 then
      D[i] = 1e-6*(3 - 5.0/3.0*(lambda_mem[i] - 3));
    else
      D[i] = 1.25e-6; //modified accordingly to literature
    end if;

    Dw[i] =D[i]*.Modelica.Math.exp(k*(1.0/Tk - 1/T_mem[i]));
    m_trans[i] = (Mv*A_mem/(0.5*t_mem)) * (Dw[i]*(C_shell[i] - C_tube[i]));
    m_trans[i] = der(C_w[i]*A_mem*t_mem);
 end for;

  C_tube = (rho_mem/M_mem) .* lambda_tube;
  C_shell = (rho_mem/M_mem) .* lambda_shell;
  T_mem = wall.Tm;

if flowConfiguration==.Modelon.ThermoFluid.Choices.FlowConfiguration.CounterFlow then
   for i in 1:n loop
     connect(shell.wall[n + 1 - i], wall.qa[i]) annotation (Line(points={{-16,50},{-16,32},{-10,32},{-10,20}}, color={191,0,0}));
     state_prim[n + 1 - i] = shell.channel.volume[i].state;
     m_trans_shell[n + 1 - i] = m_trans[i];
   end for;
else
   connect(shell.wall, wall.qa) annotation (Line(points={{-16,50},{-16,32},{-10,32},{-10,20}}, color={191,0,0}));
   state_prim = shell.channel.volume.state;
   m_trans_shell = m_trans;
end if;

  connect(wall.qb, tube.wall) annotation (Line(points={{-10,-20},{-10,-32},{0,-32},{0,-50}}, color={191,0,0}));
  connect(portB_prim, shell.portB) annotation (Line(points={{-100,40},{-66,40},{
          -66,60},{-26,60}}, color={209,60,0}));
  connect(shell.portA, portA_prim) annotation (Line(points={{10,60},{40,60},{40,
          -40},{100,-40}}, color={209,60,0}));
  connect(tube.portB, portB_sec) annotation (Line(points={{10,-60},{60,-60},{60,
          0},{100,0}}, color={209,60,0}));
  connect(portA_sec, tube.portA) annotation (Line(points={{-100,0},{-66,0},{-66,
          -60},{-26,-60}}, color={209,60,0}));
  annotation (Icon(graphics={
        Text(
          extent={{-100,-80},{100,-120}},
          textColor={0,0,127},
          textString="%name"),
        Rectangle(
          extent={{-90,60},{90,-60}},
          lineColor={0,0,0},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid),
        Rectangle(
          extent={{-90,22},{90,-22}},
          lineColor={0,0,0},
          fillPattern=FillPattern.Solid,
          fillColor={255,255,255}),
        Rectangle(
          extent={{-90,26},{90,22}},
          lineColor={0,0,0},
          fillColor={255,255,255},
          fillPattern=FillPattern.Backward),
        Rectangle(
          extent={{-90,-22},{90,-26}},
          lineColor={0,0,0},
          fillColor={255,255,255},
          fillPattern=FillPattern.Backward),
        Polygon(
          points={{-44,76},{-60,70},{-44,64},{-44,76}},
          lineColor={0,0,0},
          fillColor={0,0,0},
          fillPattern=FillPattern.Solid,
          lineThickness=0.5),
        Line(points={{60,70},{-60,70}}, color={0,0,0}),
        Line(points={{-60,-70},{60,-70}}, color={0,0,0}),
        Polygon(
          points={{42,-64},{60,-70},{42,-76},{42,-64}},
          lineColor={0,0,0},
          fillColor={0,0,0},
          fillPattern=FillPattern.Solid)}),
                                 Diagram(coordinateSystem(preserveAspectRatio=false,
          extent={{-100,-100},{100,100}})),
    Documentation(info="<html>
<p>This is a shell-and-tube type of a gas-to-gas humidifier model that includes discretized dynamic channel models which represent 
the tube and shell control volumes of the humidifier, and a dynamic wall model which acts as membrane 
through which heat and mass transfer take place. 
Co-flow and counter flow configurations are supported.</p>
<h4>Parametrization</h4>
The following humidifier membrane parameters are mapped accordingly to the control volumes (shell and tube channels) as well as the wall.
<br>
<ul>
<li><code>n</code> : number of discretizations</li>
<li><code>rho_mem</code> : membrane dry density</li>
<li><code>M_mem</code> : molar mass of dry membrane</li>
<li><code>t_mem</code> : thickness of membrane tube</li>
<li><code>D_mem</code> : inner diameter of membrane tube</li>
<li><code>D_h</code> : inner diameter of the humidifier housing</li>
<li><code>L_mem</code> : active length of membrane tube</li>
<li><code>n_tubes</code> : number of the membrane tubes</li>
<li><code>Tk</code> : material-dependent temperature constant</li>
<li><code>k</code> : material-dependent constant</li>
<li><code>kc</code> : heat transfer coefficient between fluid and substrate</li>
<li><code>A_mem</code> : pre-calculated membrane area based on the tube dimensions</li>
</ul>
<h4>Humidification</h4>
The diffusion process is driven by the water concentration gradient at the membrane boundaries.
<br>The model will first calculate the water activities of the two channels (shell and tube) 
from the gas flows, and their average as the membrane water activity according to the following formulas:
<p><img width=\"250\" src=\"modelica://FuelCell/Resources/images/formulas/humidifier1.png\"/></p>
in which p<sub>v</sub> and p<sub>sat</sub> are the partial and saturation pressure of the vapor respectively.
<br><br>The boundary water contents &lambda; for shell, tube and membrane are derived from the 
water activities according to the following:
<p><img width=\"500\" src=\"modelica://FuelCell/Resources/images/formulas/humidifier2.png\"/></p>
Water concentration on both shell and tube sides are then calculated using the provided membrane dry density and molar mass:
<p><img width=\"350\" src=\"modelica://FuelCell/Resources/images/formulas/humidifier3.png\"/></p>
The empirical constant D is determined from &lambda;<sub>mem</sub> and this constant is used for calculating 
the diffusion coefficient D<sub>w</sub> as follows:
<p><img width=\"500\" src=\"modelica://FuelCell/Resources/images/formulas/humidifier4.png\"/></p>
where k and Tk are tunable parameters which depend on the membrane material.
<br>Finally the mass transfer is given as follows:
<p><img width=\"300\" src=\"modelica://FuelCell/Resources/images/formulas/humidifier5.png\"/></p>

<h4>Assumptions</h4>
<p><ul>
<li>Gases used are ideal gases</li>
<li>Negligible potential and kinetic energy changes</li>
<li>Heat transfer only occurs across the membrane and no heat losses</li>
<li>No liquid phase species and condensation is considered</li>
<li>Diameters of control volumes equal the inner diameter</li>
</ul></p>
<h4>References</h4>
<dl><dt>Park, S.K., Choe, S.Y., and Choi, S.H.:</dt>
<p style=\"margin-left: 30px;\"><h4>Dynamic modeling and analysis of a shell-and-tube type gas-to-gas membrane humidifier for PEM fuel cell applications</h4></p>
<p style=\"margin-left: 30px;\">International Journal of Hydrogen Energy 33 (2008) 2273-2282</p>
<p style=\"margin-left: 30px;\"><a href=\"http://dx.doi.org/10.1016/j.ijhydene.2008.02.058\">DOI:10.1016/j.ijhydene.2008.02.058</a> </p>
<dl><dt>Park, S.K., Choe, S.Y., Lim, T.W., Kim, J.S., Seo, D.H., and Choi, S.H.:</dt>
<p style=\"margin-left: 30px;\"><h4>Analysis of a shell-and-tube type gas-to-gas membrane humidifier for an automotive polymer electrolyte membrane fuel cell power system</h4></p>
<p style=\"margin-left: 30px;\">International Journal of Automotive Technology, Vol. 14, No. 3, pp. 449-457 (2013)</p>
<p style=\"margin-left: 30px;\"><a href=\"http://dx.doi.org/10.1007/s12239-013-0049-4\">DOI:10.1007/s12239-013-0049-4</a> </p>
</html>", revisions="<html>
Copyright &copy; 2004-2026, MODELON AB <br /> The use of this software component is regulated by the licensing conditions for Modelon Libraries. <br />This copyright notice must, unaltered, accompany all components that are derived from, copied from, <br />or by other means have their origin from any Modelon Library.
</html>"));
end DiscretizedGasGasHumidifier;
