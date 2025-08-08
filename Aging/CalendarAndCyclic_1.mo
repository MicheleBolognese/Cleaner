within Aging;

model CalendarAndCyclic_1 "calendar and cyclic aging"
  extends .Electrification.Batteries.Core.Aging.Templates.CombinedAging(redeclare replaceable
      .Electrification.Batteries.Core.Aging.Calendar.LifetimeModel.LifetimeAgingParameter calendar, redeclare replaceable
      .Electrification.Batteries.Core.Aging.Cyclic.Continuous.CyclicContinuous cyclic,
    enable_heatport=true);
  annotation(Documentation(revisions="<html>Copyright &copy; 2004-2024, MODELON AB <br /> The use of this software component is regulated by the licensing conditions for Modelon Libraries. <br />This copyright notice must, unaltered, accompany all components that are derived from, copied from, <br />or by other means have their origin from any Modelon Library.</html>", info="<html>
<p>This is an example of a combined cyclic aging model and a calendar aging model. For information about how the contributions of these two components are combined, please refer to the <a href=\"modelica://Electrification.Batteries.Core.Aging.Templates.CombinedAging\">CombinedAging</a> template.</p>
</html>"));
end CalendarAndCyclic_1;
