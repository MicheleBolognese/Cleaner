within PEMFCModel4par.Medium;
package NASASteamH2N2O2 "NASA mixture: H2, H2O, N2, O2"

extends .FuelCell.Media.Templates.ReactionGas(
  redeclare package EquationsOfState =
      .Modelon.Media.EquationsOfState.Templates.IdealGases.NASA7 (
      data={.Modelon.Media.DataDefinitions.NasaData.H2, .Modelon.Media.DataDefinitions.NasaData.H2O, .Modelon.Media.DataDefinitions.NasaData.N2, .Modelon.Media.DataDefinitions.NasaData.O2},
      referenceChoice=.Modelon.Media.Interfaces.Types.ReferenceEnthalpy.ZeroAt25C,
      excludeEnthalpyOfFormation=false),
  redeclare package TransportProperties =
      .Modelon.Media.TransportProperties.Templates.Nasa.IdealGasNasa7TransportProperties
        (
      data={.Modelon.Media.DataDefinitions.NasaData.H2, .Modelon.Media.DataDefinitions.NasaData.H2O, .Modelon.Media.DataDefinitions.NasaData.N2, .Modelon.Media.DataDefinitions.NasaData.O2},
      substanceNames=substanceNames,
      fluidConstants=fluidConstants),
  fluidConstants={.Modelon.Media.DataDefinitions.FluidData.H2, .Modelon.Media.DataDefinitions.FluidData.H2O, .Modelon.Media.DataDefinitions.FluidData.N2, .Modelon.Media.DataDefinitions.FluidData.O2},
  gibbsFormationData={.FuelCell.Media.DataDefinitions.GibbsFormation.H2, .FuelCell.Media.DataDefinitions.GibbsFormation.H2O, .FuelCell.Media.DataDefinitions.GibbsFormation.N2, .FuelCell.Media.DataDefinitions.GibbsFormation.O2},
  reference_X=moleToMassFractions({0.25,0.25,0.25, 0.25}, MMX),
  mediumName="Gas mixture with H2O, H2, N2, O2",
  substanceNames={.Modelon.Media.DataDefinitions.NasaData.H2.name, .Modelon.Media.DataDefinitions.NasaData.H2O.name, .Modelon.Media.DataDefinitions.NasaData.N2.name, .Modelon.Media.DataDefinitions.NasaData.O2.name},
  MMX=EquationsOfState.MMX,
  R=EquationsOfState.R,
  hasCriticalData=true,
  CompressibilityType=.Modelon.Media.Interfaces.Types.Compressibility.FullyCompressible,
  analyticInverseTfromh=EquationsOfState.analyticInverseTfromh,
  ThermoStates=.Modelon.Media.Interfaces.Types.IndependentVariables.pTX,
  FixedComposition=false,
  EqReactions=0,
  nu=fill(0,0,nS),
  Hf=EquationsOfState.data.Hf,
  O_atoms={0,1,0,2},
  C_atoms={0,0,0,0},
  H_atoms={2,2,0,0});

annotation (Protection(access=Access.packageDuplicate), Documentation(info="<html>
<h4>NASASteamHydrogen</h4>
<p>NASA steam-hydrogen mixture: H2O, H2.</p>
<p>Thermodynamic properties are calculated using full 7-coefficient NASA representation (range: 200 K to 6000 K).</p>
</html>", revisions="<html>
Copyright &copy; 2004-2024, MODELON AB <br /> The use of this software component is regulated by the licensing conditions for Modelon Libraries. <br />This copyright notice must, unaltered, accompany all components that are derived from, copied from, <br />or by other means have their origin from any Modelon Library.
</html>"));
end NASASteamHydrogenNitrogen;