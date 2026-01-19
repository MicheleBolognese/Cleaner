within SAMU.HollowFiberTest2_EB_inline;

model DynamicWall "Dynamic wall with conduction heat res."
  outer .Modelon.ThermoFluid.Settings_TF settings_TF;
  extends .Modelon.ThermoFluid.Solids.Interfaces.Wall(
    qa(T(start=T0)),
    qb(T(start=T0)),
    q(T=Tm));
  extends .Modelon.ThermoFluid.Icons.WallDyn;

  parameter Integer n_channels(min=1) = 1
    "Number of parallel channels considered" annotation(Evaluate=true);

  input .Modelon.ThermoFluid.Solids.Records.BaseWallPars
                                                    props[n]
    "Geometry and material parameters, single channel, per segment"
    annotation(Dialog(enable=true));
  parameter .Modelica.Units.SI.Mass[n] m "Metal mass";
  .Modelica.Units.SI.Temperature[n] Tm(start=T0)
    "Metal mean temperature, direction along n as on side qa";
  .Modelica.Units.NonSI.Temperature_degC[n] Tm_degC=Tm - ones(n)*273.15
    "Metal mean temperature, direction along n as on side qa";
  parameter .Modelica.Units.SI.Temperature[n] T0=ones(n)*300.0
    "Metal start temperature, direction along n as on side qa";
  parameter Boolean steadyStateInit=false
    "True if derivatives should be zero at initialization";
  parameter Boolean massLessWall=false "Removes thermal inertia if true" annotation(Evaluate=true);
  parameter Boolean includeThermalResistance=true
    "Removes thermal resistance if false"                                               annotation(Evaluate=true);

  parameter Boolean TA_is_IV=true "for PbS settings only - one port should be disabled if multiple wall instances are connected in series" annotation(Dialog(tab="PbS"));
  parameter Boolean TB_is_IV=true "for PbS settings only - one port should be disabled if multiple wall instances are connected in series" annotation(Dialog(tab="PbS"));
    
  .Modelica.Units.SI.Temperature[n] TA(start=T0)=qa.T;
  .Modelica.Units.SI.Temperature[n] TB(start=T0)=qb.T;

  .Modelica.Units.SI.Power checkEnergyBalance[n] 
    "Energy balance check for each wall segment, should be zero";
initial equation
  if not massLessWall then
    if steadyStateInit then
      der(Tm) = zeros(n);
    else
      Tm = T0;
    end if;
  end if;
equation

  for i in 1:n loop
    if massLessWall then
      0.0 =  qa[i].Q_flow + qb[i].Q_flow + q[i].Q_flow
        annotation (__Modelon(ResidualEquation(iterationVariable = Tm[i])));
      checkEnergyBalance[i] = qa[i].Q_flow + qb[i].Q_flow + q[i].Q_flow;
    else
      n_channels * m[i]*props[i].Cp*der(Tm[i]) = qa[i].Q_flow + qb[i].Q_flow + q[i].Q_flow;
      checkEnergyBalance[i] = n_channels * m[i] * props[i].Cp * der(Tm[i]) - (qa[i].Q_flow + qb[i].Q_flow + q[i].Q_flow);
    end if;
    if includeThermalResistance then
      qa[i].Q_flow = n_channels * (qa[i].T - Tm[i])*2/props[i].Rw
        annotation(__Modelon(ResidualEquation(enabled = ( TA_is_IV), iterationVariable(enabled = TA_is_IV)=TA[i])));
      qb[i].Q_flow = n_channels * (qb[i].T - Tm[i])*2/props[i].Rw
        annotation(__Modelon(ResidualEquation(enabled = (TB_is_IV), iterationVariable(enabled =TB_is_IV)=TB[i])));
    else
      qa[i].T = Tm[i];
      qb[i].T = Tm[i];
    end if;
  end for;
  annotation (Documentation(revisions="<html>
Copyright &copy; 2004-2025, MODELON AB <br /> The use of this software component is regulated by the licensing conditions for Modelon Libraries. <br />This copyright notice must, unaltered, accompany all components that are derived from, copied from, <br />or by other means have their origin from any Modelon Library.
</html>"));
end DynamicWall;
