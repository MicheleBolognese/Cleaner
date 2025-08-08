within PEMFCModel.Stacks.Interfaces;
partial model AdvancedStackParameters "Advanced stack parameters"
  
  parameter Boolean enable_setting_advanced = true "If true, selectability of stack advanced parameters enabled" annotation (Dialog(tab="Advanced"));
  parameter .Modelica.Units.SI.Pressure dp_smooth = 10 "Pressure drop smoothing region around zero for anode and cathode channels" annotation (Dialog(tab="Advanced", group="Regularization",enable = enable_setting_advanced ));
  parameter .Modelica.Units.SI.MassFlowRate mflow_smooth = 1e-4 "Mass flow rate smoothing region around zero for anode and cathode channels" annotation (Dialog(tab="Advanced", group="Regularization",enable = enable_setting_advanced));

  parameter Boolean from_dp_anode = true "If true, then massflow rate is computed from pressure drop" annotation(Dialog(tab="Advanced",group="Anode side",enable = enable_setting_advanced),
        choices(choice=true "Massflow is computed from pressure",
                choice=false "Pressure drop is computed from massflow"));
  parameter Boolean positiveFlow_anode = false "Assume positive flow for upstream fluid properties" annotation(Dialog(tab="Advanced", group="Anode side",enable = enable_setting_advanced));
  parameter Boolean generateEventForReversal_anode = false "Flag for switching events for flow reversal on/off" annotation (Dialog(tab="Advanced", group="Anode side",enable = enable_setting_advanced));

  parameter Boolean from_dp_cathode = true "If true, then massflow rate is computed from pressure drop" annotation(Dialog(tab="Advanced",group="Cathode side",enable = enable_setting_advanced),
        choices(choice=true "Massflow is computed from pressure",
                choice=false "Pressure drop is computed from massflow"));
  parameter Boolean positiveFlow_cathode = false "Assume positive flow for upstream fluid properties" annotation(Dialog(tab="Advanced", group="Cathode side",enable = enable_setting_advanced));
  parameter Boolean generateEventForReversal_cathode = false "Flag for switching events for flow reversal on/off" annotation (Dialog(tab="Advanced", group="Cathode side",enable = enable_setting_advanced));

  annotation (Documentation(revisions="<html>
Copyright &copy; 2004-2024, MODELON AB <br /> The use of this software component is regulated by the licensing conditions for Modelon Libraries. <br />This copyright notice must, unaltered, accompany all components that are derived from, copied from, <br />or by other means have their origin from any Modelon Library.
</html>", info="<html>
<h4>AdvancedParameters</h4>
<p>Interface with advanced parameters.</p>
</html>"));

end AdvancedStackParameters;