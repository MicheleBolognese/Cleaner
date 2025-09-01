within SAMU.Originals;

model DistributedPipe_1_backup
  "Example with distributed pipes, different discretization numbers"
  extends .Modelon.Icons.Experiment;
  .Modelon.ThermoFluid.Sources.PressureBoundary pressureBoundary2(
    N=1,
    p=100000,
    T=298.15,
    redeclare package Medium =
        .Modelon.Media.PreDefined.IdealGases.NasaAirMixture)
              annotation (Placement(transformation(extent={{74,-86},{54,-66}})));
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
    annotation (Placement(transformation(extent={{-2,-86},{20,-66}})));

  .Modelon.ThermoFluid.Sources.MassFlowBoundary massFlowBoundary1(
    use_Th_in=true,
    m_flow=2e-3,
    T=293.15,
    redeclare package Medium =
        .Modelon.Media.PreDefined.IdealGases.NasaAirMixture)
    annotation (Placement(transformation(extent={{-50,-86},{-30,-66}})));
  .Modelica.Blocks.Sources.Ramp ramp1(
    duration=1,
    offset=298.15,
    startTime=1,
    height=20)
    annotation (Placement(transformation(extent={{-76,-42},{-56,-22}})));
  .Modelon.ThermoFluid.Sources.Environment_T environment_T1(N=10, T0=303.15)
    annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=180,
        origin={9,-48})));
equation
  connect(pipe1.portB[1], pressureBoundary2.fluidPort[1]) annotation (Line(
      points={{20,-76},{55,-76}},
      color={255,128,0},
      smooth=Smooth.None));
  connect(massFlowBoundary1.fluidPort, pipe1.portA[1]) annotation (Line(
      points={{-31,-76},{-2,-76}},
      color={255,128,0},
      smooth=Smooth.None));
  connect(ramp1.y, massFlowBoundary1.T_in)
                                         annotation (Line(
      points={{-55,-32},{-39,-32},{-39,-66}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(environment_T1.port, pipe1.q)
                                      annotation (Line(
      points={{9,-58},{9,-71}},
      color={191,0,0},
      thickness=0.5,
      smooth=Smooth.None));
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
end DistributedPipe_1_backup;
