within SAMU.Backup_preM.MembraneHumi_Mike;

partial model PartialCell_orig "General fuel cell membrane interface"

  replaceable package Medium_an =
      .FuelCell.Media.PreDefined.IdealGases.NASAReformateLong
    constrainedby .FuelCell.Media.Templates.ReactionGas "Anode Medium"      annotation(choicesAllMatching,Dialog(enable = false,tab="Predefined", group="Automatically set from component, do not override"));

  replaceable package Medium_cath =
      .FuelCell.Media.PreDefined.IdealGases.NASAMoistAir
     constrainedby .FuelCell.Media.Templates.ReactionGas "Cathode Medium"    annotation(choicesAllMatching, Dialog(enable = false,tab="Predefined", group="Automatically set from component, do not override"));

  parameter Integer N(min=1) = 1 "Number of discretization nodes" annotation(Dialog(enable = false,tab="Predefined", group="Automatically set from component, do not override"));
  constant Integer nX_an=Medium_an.nS
    "total number of mass fractions for anode";
  constant Integer nX_cath=Medium_cath.nS
    "total number of mass fractions for cathode";

  parameter .Modelica.Units.SI.Area A_cell=180e-4 "active cell area" annotation (Dialog(
      enable=false,
      tab="Predefined",
      group="Automatically set from component, do not override"));
  parameter Integer n_cell=50 "number of cells" annotation(Dialog(enable = false,tab="Predefined", group="Automatically set from component, do not override"));
  parameter .Modelica.Units.SI.HeatCapacity MCp_cell=100 "total heat capacity, PEN, substrate and electrolyte"
    annotation (Dialog(
      enable=false,
      tab="Predefined",
      group="Automatically set from component, do not override"));
  parameter .Modelica.Units.SI.CoefficientOfHeatTransfer kc=250
    "heat transfer coefficient between fluid and substrate" annotation (Dialog(
      group="Automatically set from component, do not override",
      tab="Predefined",
      enable=false));

  parameter .Modelica.Units.SI.Pressure pstart=1e5 "Pressure start value" annotation (Dialog(
      enable=false,
      tab="Predefined",
      group="Automatically set from component, do not override"));
  parameter .Modelica.Units.SI.Temperature Tstart=80 + 273 "Temperature start value" annotation (Dialog(
      enable=false,
      tab="Predefined",
      group="Automatically set from component, do not override"));

  .Modelica.Units.SI.Temperature T_cell_avg=sum(T_cell)/N "Average cell temperature";

  .Modelica.Units.SI.Temperature[N] T_cell(each start=Tstart, each fixed=true) "Substrate temperature";
  .Modelica.Units.SI.Power[N] Q_cell "heat produced in each section";
  .Modelica.Units.SI.Power Q_total "total heat produced";

  .Modelica.Units.SI.MoleFraction[N,nX_an] y_an(start=ones(N, nX_an)/nX_an) "anode gas molar fractions";
  .Modelica.Units.SI.MoleFraction[N,nX_cath] y_cath(start=ones(N, nX_cath)/nX_cath) "cathode gas molar fractions";
  // Conversion index vectors, only used in condensing case when Fluid in anode/cathode
  //   channels have different component order than membrane Medium

  parameter Boolean T_from_h = false
    "calculate anode and cathode temperature explicitly from enthalpy"                                  annotation (Dialog(
      tab="Predefined",
      group="Automatically set from component, do not override",
      enable=false));

  .Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a[N] wall(each T(nominal=500))  annotation (Placement(
        transformation(extent={{100,0},{132,32}},  rotation=0),
        iconTransformation(extent={{100,-10},{120,10}})));
  .Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_b[N] wall_an(each T(nominal=500))  annotation (Placement(
        transformation(
        origin={40,60},
        extent={{-16,-16},{16,16}},
        rotation=270), iconTransformation(
        extent={{-10,-10},{10,10}},
        rotation=270,
        origin={40,50})));
  .Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_b[N] wall_cath(each T(nominal=500))  annotation (Placement(
        transformation(
        origin={40,-60},
        extent={{-16,-16},{16,16}},
        rotation=270), iconTransformation(
        extent={{-10,-10},{10,10}},
        rotation=270,
        origin={40,-50})));
 

  .FuelCell.Interfaces.MassPort_b[N] port_an(redeclare package Medium = Medium_an)
    annotation (Placement(transformation(extent={{-60,40},{-20,80}}), iconTransformation(extent={{-40,34},{-20,54}})));
  .FuelCell.Interfaces.MassPort_b[N] port_cath(redeclare package Medium =
        Medium_cath)
    annotation (Placement(transformation(extent={{-60,-74},{-20,-34}}), iconTransformation(extent={{-40,-54},{-20,-34}})));
protected
  .Modelica.Units.SI.Temperature[N] T_an "temperature on the anode side";
  .Modelica.Units.SI.SpecificEnthalpy[N] h_an "specific enthalpy on the anode side";
  .Modelica.Units.SI.Temperature[N] T_cath "temperature on the cathode side";
  .Modelica.Units.SI.SpecificEnthalpy[N] h_cath "specific enthalpy on the cathode side";
  .Modelica.Units.SI.MassFraction[N,Medium_an.nS] X_an "sorted mass fractions on anode side";
  .Modelica.Units.SI.MassFraction[N,Medium_cath.nS] X_cath "sorted mass fractions on cathode side";

  Medium_an.ThermodynamicState[N] state_an
    "thermodynamic state on the anode side";
  Medium_cath.ThermodynamicState[N] state_cath
    "thermodynamic state on the cathode side";
  Medium_an.ReactionProperties[N] anode(p=state_an.p, T(each start=Tstart, each stateSelect=StateSelect.default)=T_an, X=X_an);
  Medium_cath.ReactionProperties[N] cathode(p=state_cath.p, T(each start=Tstart, each stateSelect=StateSelect.default)=T_cath, X=X_cath);
  Real[N] der_T_cell "Tcell derivative";

equation
  // Electrical connectors

  // Mass transfer connectors
  h_an=port_an.h;
  X_an[:,1:nX_an]=state_an.X;

  state_an=Medium_an.setState_phX(port_an.p,h_an,port_an.X);

  if Medium_an.analyticInverseTfromh or T_from_h then
      T_an = Medium_an.temperature(state_an);
  else
      h_an = Medium_an.specificEnthalpy_pTX(state_an.p,T_an,port_an.X);
  end if;

  h_cath=port_cath.h;
  X_cath[:,1:nX_cath]=state_cath.X;

  state_cath=Medium_cath.setState_phX(port_cath.p,h_cath,port_cath.X);

  if Medium_cath.analyticInverseTfromh or T_from_h then
    T_cath = Medium_cath.temperature(state_cath);
  else
    h_cath = Medium_cath.specificEnthalpy_pTX(port_cath.p,T_cath,port_cath.X);
  end if;

  // Heat transfer connectors, assume heat transfer area=2*cell area
  wall.T = T_cell;
  wall_an.Q_flow = kc*n_cell*A_cell/(N)*(wall_an.T - T_cell);
  wall_cath.Q_flow = kc*n_cell*A_cell/(N)*(wall_cath.T - T_cell);

  der_T_cell=(wall_an.Q_flow + wall_cath.Q_flow + wall.Q_flow + Q_cell)*N/MCp_cell;
  MCp_cell/(N)*der(T_cell) = wall_an.Q_flow + wall_cath.Q_flow + wall.Q_flow + Q_cell;

  for i in 1:N loop
    y_an[i,:]  = anode[i].Z;
    y_cath[i,:]= cathode[i].Z;
  end for;


  annotation (Icon(
      coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,100}}),
        graphics={Rectangle(
          extent={{-100,40},{100,-40}},
          lineColor={0,0,255},
          fillColor={135,135,135},
          fillPattern=FillPattern.CrossDiag), Text(
          extent={{-100,8},{100,-12}},
          lineColor={255,255,255},
          textString="%name")}),
    Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,
            100}}),     graphics),
    Documentation(info="<html>
<h4>PartialCell</h4>
<p>Partial cell membrane base model with the connectors listed below and heat transfer equations. To create a complete model, equations for the electrochemical properties must be added.</p>
<ul>
<li>Electrical pins (positive and negative)</li>
<li>Heat port connector (wall), discretized in the flow direction to describe the temperature profile in the flow direction of the membrane </li>
<li>Heat port connectors for heat transfer between membrane and flow channel</li>
<li>Component mass transfer connectors of reaction gas type (ideal gas) to the anode and cathode sides</li>
</ul>
<h4>Parametrization</h4>
<p>All the parameters of the model have predefined values which the user should not change. The discretization of the model is defined in the following way:</p>
<ul>
<li>The Parameter <span style=\"font-family: Courier New;\">N</span> describes the number of discretizations in the flow directions (membrane divided in N segments) </li>
<li>The parameter <span style=\"font-family: Courier New;\">n_cell</span> describes the number of cells in one membrane; the membrane can be used to describe a single cell membrane, or a number of cells (for lumped stacks)</li>
</ul>
</html>", revisions="<html>
Copyright &copy; 2004-2025, MODELON AB <br /> The use of this software component is regulated by the licensing conditions for Modelon Libraries. <br />This copyright notice must, unaltered, accompany all components that are derived from, copied from, <br />or by other means have their origin from any Modelon Library.
</html>"));
end PartialCell_orig;
