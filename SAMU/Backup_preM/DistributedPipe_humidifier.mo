within SAMU.Backup_preM;

model DistributedPipe_humidifier
  "Example with distributed pipes, different discretization numbers"
  //extends .Modelon.Icons.Experiment;
  //extends .FuelCell.Icons.Pipe;  
  .Modelon.ThermoFluid.FlowChannels.DistributedPipe pipe1(
    Dhyd=1e-2,
    L=2,
    redeclare model HeatTransfer =
        .Modelon.ThermoFluid.FlowChannels.HeatTransfer.SinglePhase.ConstantCoefficient
        (                                                            alpha=100),
    n=10,
    redeclare model Friction =
        .Modelon.ThermoFluid.FlowChannels.PipeResistances.SinglePhase.DarcyWeisbach
        (zeta_A=1),
    redeclare package Medium =
        .Modelon.Media.PreDefined.IdealGases.NasaAirMixture,
    D=1e-2,
    A=2e-4,
    A_heat=0.1,
    useHeatTransfer=true,
    NA=1,
    NB=1)
    annotation (Placement(transformation(extent={{-2.0,-86.0},{20.0,-66.0}},rotation = 0.0,origin = {0.0,0.0})));
  .Modelon.ThermoFluid.Sources.Environment_T environment_T1(N=10, T0=303.15)
    annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=180,
        origin={9,-48})));
    .Modelon.ThermoFluid.Interfaces.VolumePort volumePort annotation(Placement(transformation(extent = {{-54.0,-86.0},{-34.0,-66.0}},origin = {0.0,0.0},rotation = 0.0)));
    .Modelon.ThermoFluid.Interfaces.FlowPort flowPort annotation(Placement(transformation(extent = {{58.0,-86.0},{78.0,-66.0}},origin = {0.0,0.0},rotation = 0.0)));
equation
  connect(environment_T1.port, pipe1.q)
                                      annotation (Line(
      points={{9,-58},{9,-71}},
      color={191,0,0},
      thickness=0.5,
      smooth=Smooth.None));
    connect(pipe1.portA[1],volumePort) annotation(Line(points = {{-2,-76},{-44,-76}},color = {255,128,0}));
    connect(pipe1.portB[1],flowPort) annotation(Line(points = {{20,-76},{68,-76}},color = {255,128,0}));
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
