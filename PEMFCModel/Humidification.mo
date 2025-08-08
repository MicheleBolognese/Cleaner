within PEMFCModel;
model Humidification
    
  replaceable package Medium = .PEMFCModel.Medium.NASALowQualityHydrogen constrainedby .FuelCell.Media.Templates.ReactionGas "Medium";
  package WaterMedium = .FuelCell.Media.PreDefined.TwoPhase.WaterIF97 "Water medium";
        
  constant Modelica.Units.SI.Pressure p_ref = 1.01325e5 "Standard pressure";
  constant Modelica.Units.SI.Temperature T_ref = 273.15 "Standard temperature";
  
  .Modelica.Blocks.Interfaces.RealInput p_in(unit = "Pa") "Gas inlet pressure" annotation (Placement(
          transformation(extent={{-140.0,40.0},{-100.0,80.0}},rotation = 0.0,origin = {0.0,0.0})));
  .Modelica.Blocks.Interfaces.RealInput V_flow_dry_in(unit = "NLPM") "Dry gas inlet volumetric flow rate" annotation (Placement(
          transformation(extent={{-140.0,-80.0},{-100.0,-40.0}},rotation = 0.0,origin = {0.0,0.0})));
  .Modelica.Blocks.Interfaces.RealOutput m_flow_wet_in(unit = "kg/s") = m_flow_dry_in/(1-x_wet_in[i_H2O]) "Wet gas inlet mass flow rate" annotation (Placement(
          transformation(extent={{100.25,21.0},{140.25,61.0}},rotation = 0.0,origin = {0.0,0.0})));
  .Modelica.Blocks.Interfaces.RealOutput x_wet_in[Medium.nS]= Medium.moleToMassFractions(y_wet_in, Medium.MMX) "Wet gas inlet mass-based composition" annotation (Placement(
          transformation(extent={{103.75,-53.0},{143.75,-13.0}},rotation = 0.0,origin = {0.0,0.0})));
        
  parameter Modelica.Units.SI.Temperature T_in = T_ref "Gas inlet temperature";
  parameter Modelica.Units.SI.MoleFraction y_dry_in[Medium.nS] = Medium.massToMoleFractions(Medium.reference_X,Medium.MMX) "Dry gas inlet molar-based composition";
  parameter Modelica.Units.SI.Temperature T_dew = 52.5 + 273.13 "Wet gas dew point";      
  
    
  Modelica.Units.SI.MolarMass MM_dry = y_dry_in*Medium.MMX "Dry gas molar mass";
  Modelica.Units.SI.Density rho_dry = p_ref*MM_dry/.Modelica.Constants.R/T_ref "Dry gas density in standard conditions";
  Modelica.Units.SI.MassFlowRate m_flow_dry_in = V_flow_dry_in*1e-3/60*rho_dry "Dry gas inlet mass flow rate";
  
  Modelica.Units.SI.MoleFraction y_H2O = WaterMedium.saturationPressure_TX(T_dew)/p_in "H2O molar fraction in wet gas";
  Modelica.Units.SI.DimensionlessRatio RH = y_H2O*p_in/WaterMedium.saturationPressure_TX(T_in);
  Modelica.Units.SI.MoleFraction y_wet_in[Medium.nS] = {if i == i_H2O then y_H2O else y_dry_in[i]*(1-y_H2O) for i in 1:Medium.nS};
  
  protected
  
    parameter Integer i_H2O = Medium.substanceIndex("H2O") "Index of H2O in the medium";
  
  initial equation 
    
     assert(i_H2O > 0, "H2O must be present in the medium", level = AssertionLevel.error);
  
   annotation(Icon(coordinateSystem(preserveAspectRatio = false,extent = {{-100.0,-100.0},{100.0,100.0}}),graphics = {Rectangle(lineColor={0,0,0},fillColor={230,230,230},fillPattern=FillPattern.Solid,extent={{-100.0,-100.0},{100.0,100.0}}),Text(lineColor={0,0,255},extent={{-150,150},{150,110}},textString="%name")}));
    
end Humidification;