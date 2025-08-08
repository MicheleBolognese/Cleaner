within PEMFCModel.Medium;
package NASAMixture "NASA industrial quality hydrogen: H2, CH4, NH3, CO, CO2, H2O, N2, O2"

extends .FuelCell.Media.Templates.ReactionGas(
  redeclare package EquationsOfState =
    .Modelon.Media.EquationsOfState.Templates.IdealGases.NASA7 (
      data = {.Modelon.Media.DataDefinitions.NasaData.H2,.Modelon.Media.DataDefinitions.NasaData.CH4,.Modelon.Media.DataDefinitions.NasaData.NH3,.Modelon.Media.DataDefinitions.NasaData.CO,.Modelon.Media.DataDefinitions.NasaData.CO2,.Modelon.Media.DataDefinitions.NasaData.H2O,.Modelon.Media.DataDefinitions.NasaData.N2,.Modelon.Media.DataDefinitions.NasaData.O2},
      referenceChoice=.Modelon.Media.Interfaces.Types.ReferenceEnthalpy.ZeroAt25C,
      excludeEnthalpyOfFormation=false),
  redeclare package TransportProperties =
    .Modelon.Media.TransportProperties.Templates.Nasa.IdealGasNasa7TransportProperties
      (data = {.Modelon.Media.DataDefinitions.NasaData.H2,.Modelon.Media.DataDefinitions.NasaData.CH4,.Modelon.Media.DataDefinitions.NasaData.NH3,.Modelon.Media.DataDefinitions.NasaData.CO,.Modelon.Media.DataDefinitions.NasaData.CO2,.Modelon.Media.DataDefinitions.NasaData.H2O,.Modelon.Media.DataDefinitions.NasaData.N2,.Modelon.Media.DataDefinitions.NasaData.O2},
      substanceNames=substanceNames,
      fluidConstants= fluidConstants),
  fluidConstants = {.Modelon.Media.DataDefinitions.FluidData.H2,.Modelon.Media.DataDefinitions.FluidData.CH4,.Modelon.Media.DataDefinitions.FluidData.NH3,.Modelon.Media.DataDefinitions.FluidData.CO,.Modelon.Media.DataDefinitions.FluidData.CO2,.Modelon.Media.DataDefinitions.FluidData.H2O,.Modelon.Media.DataDefinitions.FluidData.N2,.Modelon.Media.DataDefinitions.FluidData.O2},
  gibbsFormationData={.PEMFCModel.Medium.GibbsFormationNH3.H2, .PEMFCModel.Medium.GibbsFormationNH3.CH4, .PEMFCModel.Medium.GibbsFormationNH3.NH3, .PEMFCModel.Medium.GibbsFormationNH3.CO, .PEMFCModel.Medium.GibbsFormationNH3.CO2, .PEMFCModel.Medium.GibbsFormationNH3.H2O, .PEMFCModel.Medium.GibbsFormationNH3.N2, .PEMFCModel.Medium.GibbsFormationNH3.O2},
  reference_X=moleToMassFractions({0.2766,0.1,0.1196,0.0180,0.0733,0.1953,0.2055,0.0117},MMX),
  mediumName="Industrial quality H2 with CH4, NH3, CO, CO2, H2O, N2, O2",
  substanceNames={.Modelon.Media.DataDefinitions.NasaData.H2.name,.Modelon.Media.DataDefinitions.NasaData.CH4.name,.Modelon.Media.DataDefinitions.NasaData.NH3.name,.Modelon.Media.DataDefinitions.NasaData.CO.name,.Modelon.Media.DataDefinitions.NasaData.CO2.name,.Modelon.Media.DataDefinitions.NasaData.H2O.name,.Modelon.Media.DataDefinitions.NasaData.N2.name,.Modelon.Media.DataDefinitions.NasaData.O2.name},
  MMX=EquationsOfState.MMX,
  R=EquationsOfState.R,
  hasCriticalData=true,
  CompressibilityType = .Modelon.Media.Interfaces.Types.Compressibility.FullyCompressible,
  analyticInverseTfromh = EquationsOfState.analyticInverseTfromh,
  ThermoStates = .Modelon.Media.Interfaces.Types.IndependentVariables.pTX,
  FixedComposition=false,
  EqReactions=0,
  nu=fill(0,0,nS),
  Hf=EquationsOfState.data.Hf,
  O_atoms={0,0,0,1,2,1,0,2},
  C_atoms={0,1,0,1,1,0,0,0},
  H_atoms={2,4,3,0,0,2,0,0});

  constant Real S_ionic[:]={1,0,0,0,-1,0,0}
  "Stoich vector for anode side electrochemistry";
  constant Real S_pox[:]={2,-1,1,0,0,0,-0.5}
  "Stoich vector for partial oxidation";
  constant Real S_comb[:,:]=
    [0,0,0,0,0,0,0; 0,1,1,1,0,0,0; 1,2,0,0,1,0,0; 0,0,0,0,0,1,0; -0.5,-2,-0.5,0,0,0,1]
  "Stoich matrix for combustion into moist air";
  constant Real C_vect[:] = {0,1,0,0,0,0,0}
  "Carbon number for calculating St2C-ratio";

annotation (Protection(access=Access.packageDuplicate), Documentation(info="<html>
<h4>NASAMixture</h4>
<p>NASA hydrogen rich gas: H2, CH4, NH3, CO, CO2, H2O, N2, O2.</p>
<p>Thermodynamic properties are calculated using full 7-coefficient NASA representation (range: 200 K to 6000 K).</p>
</html>", revisions="<html>
Copyright &copy; 2004-2024, MODELON AB <br /> The use of this software component is regulated by the licensing conditions for Modelon Libraries. <br />This copyright notice must, unaltered, accompany all components that are derived from, copied from, <br />or by other means have their origin from any Modelon Library.
</html>"));
end NASAMixture;
