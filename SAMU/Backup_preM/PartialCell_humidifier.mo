within SAMU.Backup_preM;

partial model PartialCell_humidifier "General fuel cell membrane interface"

  replaceable package Medium_an = .FuelCell.Media.PreDefined.IdealGases.NASAReformate constrainedby 
    .FuelCell.Media.Templates.ReactionGas "M;edium on the dry side" annotation(choicesAllMatching,Dialog(enable = enable_setting));

  replaceable package Medium_cath = .FuelCell.Media.PreDefined.IdealGases.NASAMoistAir constrainedby
    .FuelCell.Media.Templates.ReactionGas "Medium on the wet side" annotation(choicesAllMatching,Dialog(enable = enable_setting));
   
  parameter Boolean enable_setting = true " = true to set cell parameters";  

  parameter Integer N(min = 1) = 1 "Number of along-the-channel discretization nodes" annotation(Dialog(enable = enable_setting));
  parameter Integer n_cell(min = 1) = 1 "Number of series-connected cell" annotation(Dialog(enable = enable_setting));
  
  parameter .Modelica.Units.SI.Area A_cell = 180e-4 "Cell active area" annotation(Dialog(enable = enable_setting));
  parameter .Modelica.Units.SI.HeatCapacity C_cell = 100 "Cell heat capacity, including substrate and electrolyte" annotation(Dialog(enable = enable_setting));
  input .Modelica.Units.SI.CoefficientOfHeatTransfer h_conv_an[N] "Heat transfer coefficient between fluid and anode substrate" annotation(Dialog(enable = enable_setting));
  input .Modelica.Units.SI.CoefficientOfHeatTransfer h_conv_cath[N] "Heat transfer coefficient between fluid and cathode substrate" annotation(Dialog(enable = enable_setting));

  parameter .Modelica.Units.SI.Pressure pstart = 1e5 "Pressure initialization value" annotation(Dialog(enable = enable_setting));
  parameter .Modelica.Units.SI.Temperature Tstart = 80 + 273.15 "Temperature initialization value" annotation(Dialog(enable = enable_setting));
  
  parameter Boolean T_from_h = false "Calculate anode and cathode temperature explicitly from enthalpy" annotation(Dialog(enable = enable_setting)); 
  
  .Modelica.Units.SI.Temperature T_cell[N](each start = Tstart, each fixed = true) "Cell substrate temperature in each node";
  .Modelon.Media.Units.DerTemperatureByTime dTdt[N] "Time-derivative of cell substrate temperature in each node";
  final .Modelica.Units.SI.Temperature T_cell_avg = sum(T_cell)/N "Average cell substrate temperature";
//   .Modelica.Units.SI.Current I_stack = pin_n.i "Stack current";
//   .Modelica.Units.SI.CurrentDensity j(min = 0) = I_stack/A_cell "Current density referred to cell active area";
//   .Modelica.Units.SI.CurrentDensity j_ionic[N](each min = 0, each start = 1) "Current density applied to each node";
//   .Modelica.Units.SI.Voltage V_cell "Single cell voltage";
//   .Modelica.Units.SI.Voltage V_stack "Stack voltage";
//   .Modelica.Units.SI.Power P_stack "Stack electrical power output";
  
  .Modelica.Units.SI.HeatFlowRate Q_stack "Stack heat production, coming from enthalpy flows and generated power";
  .Modelica.Units.SI.HeatFlowRate Q_wall_an_stack "Stack heat flow rate between anode channel and cell substrate";
  .Modelica.Units.SI.HeatFlowRate Q_wall_cath_stack "Stack heat flow rate between cathode channel and cell substrate";
  .Modelica.Units.SI.HeatFlowRate Q_wall_stack "Stack heat flow rate between cell substrate and external object";
  .Modelica.Units.SI.MoleFraction y_an[N, nS_an](start = ones(N, nS_an)/nS_an) "Anodic gas molar-based composition in each node";
  .Modelica.Units.SI.MoleFraction y_cath[N, nS_cath](start = ones(N, nS_cath)/nS_cath) "Cathodic gas molar-based composition in each node";

  .Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a wall[N](each T(nominal = 500)) annotation (Placement(
   transformation(extent={{100,0},{132,32}}, rotation=0),     iconTransformation(extent={{100,-10},{120,10}})));
  .Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_b wall_an[N](each T(nominal = 500)) annotation (Placement(
        transformation(
        origin={40,60},
        extent={{-16,-16},{16,16}},
        rotation=270), iconTransformation(
        extent={{-10,-10},{10,10}},
        rotation=270,
        origin={40,50})));
  .Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_b wall_cath[N](each T(nominal = 500)) annotation (Placement(
        transformation(
        origin={40,-60},
        extent={{-16,-16},{16,16}},
        rotation=270), iconTransformation(
        extent={{-10,-10},{10,10}},
        rotation=270,
        origin={40,-50})));
//   .Modelica.Electrical.Analog.Interfaces.PositivePin pin_p "Stack positive terminal" annotation (Placement(transformation(extent={{-100.0,-60.0},{-80.0,-40.0}},
//           rotation=0.0,origin = {0.0,0.0})));
//   .Modelica.Electrical.Analog.Interfaces.NegativePin pin_n "Stack negative terminal" annotation (Placement(transformation(extent={{-99.68694318907956,40.0},{-79.68694318907956,60.0}}, rotation=
//             0.0,origin = {0.0,0.0})));
  .FuelCell.Interfaces.MassPort_b port_an[N](redeclare package Medium = Medium_an) annotation (Placement(transformation(extent={{-60,40},{-20,80}}), iconTransformation(extent={{-40,34},{-20,54}}))); 
  .FuelCell.Interfaces.MassPort_b port_cath[N](redeclare package Medium = Medium_cath) annotation (Placement(transformation(extent={{-60,-74},{-20,-34}}), iconTransformation(extent={{-40,-54},{-20,-34}})));
   
   final parameter Integer nS_an = Medium_an.nS "Number of chemical species in the anodic medium" annotation(Dialog(enable = false));
   final parameter Integer nS_cath = Medium_cath.nS "Number of chemical species in the cathodic medium" annotation(Dialog(enable = false));
//   .Modelica.Units.SI.Power P_cell[N] "Single cell electrical power output in each node";
//  .Modelica.Units.SI.HeatFlowRate Q_cell[N] "Sigle cell heat production in each node, coming from enthalpy flows and generated power";
  .Modelica.Units.SI.HeatFlowRate Q_wall_an_cell[N] "Single cell heat flow rate between anode channel and cell substrate in each node";
  .Modelica.Units.SI.HeatFlowRate Q_wall_cath_cell[N] "Single cell heat flow rate between cathode channel and cell substrate in each node";
  .Modelica.Units.SI.HeatFlowRate Q_wall_cell[N] "Single cell heat flow rate between cell substrate and external object in each node";
  .Modelica.Units.SI.Temperature T_an[N] "Dry side gas temperature in each node at channel-membrane interface";
  .Modelica.Units.SI.Pressure p_an[N](each start = pstart) "Dry side gas pressure in each node at channel-membrane interface";
  .Modelica.Units.SI.SpecificEnthalpy h_an[N] "Dry side gas specific enthalpy in each node at channel-membrane interface";
  .Modelica.Units.SI.MassFraction X_an[N, nS_an] "Dry side gas mass fractions in each node at channel-membrane interface";
  .Modelica.Units.SI.Temperature T_cath[N] "Wet side gas temperature in each node at channel-membrane interface";
  .Modelica.Units.SI.Pressure p_cath[N](each start = pstart) "Wet side gas pressure in each node at channel-membrane interface";
  .Modelica.Units.SI.SpecificEnthalpy h_cath[N] "Wet side gas specific enthalpy in each node at channel-membrane interface";
  .Modelica.Units.SI.MassFraction X_cath[N, nS_cath] "Wet side gas mass fractions in each node at channel-membrane interface";
  Medium_an.ReactionProperties anode[N](p = p_an, T(each start = Tstart, each stateSelect = StateSelect.default) = T_an, X = X_an);
  Medium_cath.ReactionProperties cathode[N](p = p_cath, T(each start = Tstart, each stateSelect = StateSelect.default) = T_cath, X = X_cath);

protected

  Medium_an.ThermodynamicState state_an[N] "Dry side gas thermodynamic state in each node at channel-membrane interface";
  Medium_cath.ThermodynamicState state_cath[N] "Wet side gas thermodynamic state in each node at channel-membrane interface"; 

equation
//   // Electrical connectors
//   pin_p.v - pin_n.v = V_stack;

  // Mass transfer connectors
  p_an = port_an.p;
  h_an = port_an.h;
  X_an = port_an.X;

  state_an = Medium_an.setState_phX(p_an,h_an,X_an);

  if Medium_an.analyticInverseTfromh or T_from_h then
      
    T_an = Medium_an.temperature(state_an);
  
  else
  
    h_an = Medium_an.specificEnthalpy_pTX(p_an,T_an,X_an);
  
  end if;

  p_cath = port_cath.p;
  h_cath = port_cath.h;
  X_cath = port_cath.X;

  state_cath = Medium_cath.setState_phX(p_cath,h_cath,X_cath);

  if Medium_cath.analyticInverseTfromh or T_from_h then
    
    T_cath = Medium_cath.temperature(state_cath);
  
  else
    
    h_cath = Medium_cath.specificEnthalpy_pTX(p_cath,T_cath,X_cath);
  
  end if;

  // Heat transfer connectors, assuming that per each side the heat transfer area is equal to the cell active area
  wall.T = T_cell;
  wall_an.Q_flow = n_cell*Q_wall_an_cell;
  wall_cath.Q_flow = n_cell*Q_wall_cath_cell;
  wall.Q_flow = n_cell*Q_wall_cell;
  Q_wall_an_cell = h_conv_an.*(wall_an.T - T_cell)*A_cell/N;
  Q_wall_cath_cell = h_conv_cath.*(wall_cath.T - T_cell)*A_cell/N;
  
  //Energy balance equation on the single cell per each cell node 
  C_cell/N*der(T_cell) = Q_wall_an_cell + Q_wall_cath_cell + Q_wall_cell /*+ Q_cell*/;
  dTdt = der(T_cell);
  
  for i in 1:N loop
    
    y_an[i, :] = Medium_an.massToMoleFractions(X_an[i,:], Medium_an.MMX);
    y_cath[i, :] = Medium_cath.massToMoleFractions(X_cath[i,:], Medium_cath.MMX);
  
  end for;

//   //Single cell quantities
//   V_cell = V_stack/n_cell;
  
//   //Stack quantities
//   P_stack = V_stack*I_stack;
  Q_wall_an_stack = sum(wall_an.Q_flow);
  Q_wall_cath_stack = sum(wall_cath.Q_flow);
  Q_wall_stack = sum(wall.Q_flow);

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
Copyright &copy; 2004-2024, MODELON AB <br /> The use of this software component is regulated by the licensing conditions for Modelon Libraries. <br />This copyright notice must, unaltered, accompany all components that are derived from, copied from, <br />or by other means have their origin from any Modelon Library.
</html>"));

end PartialCell_humidifier;
