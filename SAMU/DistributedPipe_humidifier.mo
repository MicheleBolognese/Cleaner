within SAMU;

model DistributedPipe_humidifier
  "Example with distributed pipes, different discretization numbers"
    .Modelon.ThermoFluid.Interfaces.VolumePort volumePort annotation(Placement(transformation(extent = {{-54.0,-86.0},{-34.0,-66.0}},origin = {0.0,0.0},rotation = 0.0)));
    .Modelon.ThermoFluid.Interfaces.FlowPort flowPort annotation(Placement(transformation(extent = {{58.0,-86.0},{78.0,-66.0}},origin = {0.0,0.0},rotation = 0.0)));
    .FuelCell.Pipes.FlowChannel_mass flowChannel_mass annotation(Placement(transformation(extent = {{11.31,-10.0},{-11.31,10.0}},origin = {10.0,-76.0},rotation = 180.0)));
    .FuelCell.Pipes.FlowChannel_mass flowChannel_mass3 annotation(Placement(transformation(extent = {{-11.2,10.84},{11.2,-10.84}},origin = {10.66,-130.0},rotation = -180.0)));
    .FuelCell.Interfaces.WaterVolumePort waterVolumePort annotation(Placement(transformation(extent = {{60.66,-137.34},{75.34,-122.66}},origin = {0.0,0.0},rotation = 0.0)));
    .FuelCell.Interfaces.WaterFlowPort waterFlowPort annotation(Placement(transformation(extent = {{-51.34,-137.34},{-36.66,-122.66}},origin = {0.0,0.0},rotation = 0.0)));
    .SAMU.Empirical_humidificatore_def empirical_humidificatore_def annotation(Placement(transformation(extent = {{2.0,-116.0},{22.0,-96.0}},origin = {0.0,0.0},rotation = 0.0)));
equation
    connect(flowChannel_mass.portA,volumePort) annotation(Line(points = {{-0.18,-76},{-44,-76}},color = {255,128,0}));
    connect(flowChannel_mass.portB,flowPort) annotation(Line(points = {{20.18,-76},{68,-76}},color = {255,128,0}));
    connect(flowChannel_mass3.portA,waterVolumePort) annotation(Line(points = {{20.74,-130},{68,-130}},color = {0,0,255}));
    connect(flowChannel_mass3.portB,waterFlowPort) annotation(Line(points = {{0.58,-130},{-44,-130}},color = {0,0,255}));
    connect(flowChannel_mass.port,empirical_humidificatore_def.port_an) annotation(Line(points = {{6.61,-80.4},{6.61,-91},{9,-91},{9,-101.6}},color = {0,0,255}));
    connect(empirical_humidificatore_def.port_cath,flowChannel_mass3.port) annotation(Line(points = {{9,-110.4},{9,-117.81},{14.02,-117.81},{14.02,-125.23}},color = {0,0,255}));
    connect(flowChannel_mass.wall,empirical_humidificatore_def.wall_an) annotation(Line(points = {{14.52,-81},{14.52,-91},{16,-91},{16,-101}},color = {191,0,0}));
    connect(empirical_humidificatore_def.wall_cath,flowChannel_mass3.wall) annotation(Line(points = {{16,-111},{16,-117.79},{6.18,-117.79},{6.18,-124.58}},color = {191,0,0}));
  annotation (
    Documentation(info="<html>
<p><h4>Purpose</h4></p>
<p>Test the DistributedPipe component by observing temperature behaviour along the pipe when inlet temperature is ramped and inlet mass flow rate is constant.</p>
<p><h4>Test criteria</h4></p>
<p>Simulates for 5 seconds with a match for the trajectory of the temperatures specified below, where the latter is delayed due to medium residence time in the pipe.</p>
<p><h4>Details</h4></p>
<p>Compare trajectories for:</p>
<p><ul>
<li>pipe.volume[1].T</li>
<li>pipe.volume[2].T</li>
<li>pipe.volume[3].T</li>
<li>pipe.volume[4].T</li>
</ul></p>
</html>"),
    experiment(StopTime=5, Tolerance=1e-05));
end DistributedPipe_humidifier;
