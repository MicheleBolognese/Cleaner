within SAMU;

partial model PartialMembrane_hgumi "Cell membrane with underlying membrane equations and loss models added"
  extends SAMU.PartialCellTransport_humidifier;

  parameter String an_names[:]={"H2"} "Required anode species | Used for index references (iAn)";
  parameter String cath_names[:]={"O2","H2O"} "Required anode species | Used for index references (iCath)";



  .Modelica.Units.SI.Power sum_Q_cell "Heat of Cell";




 



  final parameter Integer[size(an_names, 1)] iAn=Medium_an.substanceIndexVector(an_names) "Indices of required anode species";
  final parameter Integer[size(cath_names, 1)] iCath=Medium_cath.substanceIndexVector(cath_names) "Indices of required cathode species";

initial equation
  assert(
    min(iAn) > 0,
    "Membrane: The selected media for the anode side does not contain the required substances: see 'an_names'",
    level=AssertionLevel.error);
  assert(
    min(iCath) > 0,
    "Membrane: The selected media for the anode side does not contain the required substances: see 'cath_names'",
    level=AssertionLevel.error);

equation


  // Reaction flow terms to mass transfer ports
  for i in 1:N loop
    for j in 1:nX_an loop
      port_an[i].mX_flow[j] =  mX_flow_an_transport[i, j];
    end for;

    port_an[i].m_flow = sum(mX_flow_an_transport[i, :]);
    port_an[i].H_flow =  sum(anode[i].h_component .* mX_flow_an_transport[i,
      :]);

    for j in 1:nX_cath loop
      port_cath[i].mX_flow[j] =  mX_flow_cath_transport[i, j];
    end for;

    port_cath[i].m_flow =  sum(mX_flow_cath_transport[i, :]);
    port_cath[i].H_flow =  sum(cathode[i].h_component .*
      mX_flow_cath_transport[i, :]);
  end for;

  // Summary calculations
  // Difference between Q_total and Q_cell is due to inconsistency in CellReaction function but should be close and can be used as a checksum
  Q_cell = port_an.H_flow + port_cath.H_flow;
  sum_Q_cell = sum(Q_cell);
  Q_total = sum(wall_an.Q_flow) - sum(wall_cath.Q_flow) - sum(wall.Q_flow);

 

  annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(coordinateSystem(preserveAspectRatio=false)),
    Documentation(info="<html>
<h4>PartialMembrane</h4>
<p>Partial cell membrane model with mass and energy balances defined based on variable selections.</p>
<p>To create a complete model:</p>
<ul>
<li>Proper replaceable models should be selected and all parameters properly tuned (e.g., for PEMFC vs SOFC).</li>
<li>One additional equation is required to define the open circuit voltage (e.g, nernst equation - E0) in the extending models</li>
</ul>
</html>", revisions="<html>
Copyright &copy; 2004-2025, MODELON AB <br /> The use of this software component is regulated by the licensing conditions for Modelon Libraries. <br />This copyright notice must, unaltered, accompany all components that are derived from, copied from, <br />or by other means have their origin from any Modelon Library.
</html>"));
end PartialMembrane_hgumi;
