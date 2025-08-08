within PEMFCModel.Stacks.Interfaces;
partial model AdvancedCoolingParameters "Set of advanced parameters for cooling channel models"
  
  parameter Boolean enable_setting_cooling_advanced = true "If true, selectability of cooling chanell advanced parameters enabled"  annotation(Dialog(tab="Advanced"));
  parameter Boolean from_dp_cooling = true "If true, then massflow rate is computed from pressure"  annotation(Dialog(tab="Advanced",group="Cooling channel",enable = enable_setting_cooling_advanced), choices(choice=true "Massflow is computed from pressure",
                choice=false "Pressure drop is computed from massflow"));  
  parameter Boolean positiveFlow_cooling = false  "Assume positive flow for upstream properties" annotation(Dialog(tab="Advanced",group="Cooling channel",enable = enable_setting_cooling_advanced));
  parameter Boolean generateEventForReversal_cooling = false "Flag for switching events for flow reversal on/off" annotation (Dialog(tab="Advanced", group="Cooling channel",enable = enable_setting_cooling_advanced));
  parameter .Modelica.Units.SI.Pressure dp_smooth_cooling = 10 "Pressure drop smoothing region around zero for cooling channel" annotation (Dialog(tab="Advanced", group="Regularization",enable = enable_setting_cooling_advanced));
  parameter .Modelica.Units.SI.MassFlowRate mflow_smooth_cooling = 0.001 "Massflow smoothing region around zero for cooling channel" annotation (Dialog(tab="Advanced", group="Regularization",enable = enable_setting_cooling_advanced));
 
  annotation (Documentation(revisions="<html>
Copyright &copy; 2004-2024, MODELON AB <br /> The use of this software component is regulated by the licensing conditions for Modelon Libraries. <br />This copyright notice must, unaltered, accompany all components that are derived from, copied from, <br />or by other means have their origin from any Modelon Library.
</html>", info="<html>
<h4>AdvancedCoolingParameters</h4>
<p>Interface with advanced cooling parameters.</p>
</html>"));

end AdvancedCoolingParameters;