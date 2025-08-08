within PEMFCModel.Membrane.Losses.Templates;
partial model BaseContaminantsLoss "Template for voltage loss due to contaminants"
  
  extends PEMFCModel.Membrane.Losses.Templates.BaseLoss;
 
  replaceable package Medium_an = .FuelCell.Media.PreDefined.IdealGases.NASAReformate constrainedby FuelCell.Media.Templates.ReactionGas annotation (Dialog(enable = enableInternal, tab = "Internal", group = "Input parameters"));
  
  parameter String contaminantsSpecies[:] = fill("", 0) "Name of contaminant species affecting voltage loss, ex: {'CO','NH3'}" annotation (Dialog(group = "Input parameters"));
  
  // Inputs to be determined! (Which parameters mostly affect cell performance degradation due to contaminants?)

  input .Modelica.Units.SI.Temperature T_cell[N] "Cell substrate temperature in each node" annotation (Dialog(enable = enableInternal, tab = "Internal", group = "Input variables")); 
  input .Modelica.Units.SI.Pressure p_an_partial[N, Medium_an.nS] "Partial pressure of each chemical species in each node at anode-membrane interface" annotation (Dialog(enable = enableInternal, tab = "Internal", group = "Input variables"));
  //input .Modelica.Units.SI.Pressure[N, nCath] p_cath_partial;
  input .Modelica.Units.SI.MoleFraction y_an[N, Medium_an.nS] "Molar-based composition at anode-membrane interface" annotation (Dialog(enable = enableInternal, tab = "Internal", group = "Input variables"));
  
  final parameter Integer nCont = size(contaminantsSpecies, 1) "Number of considered contaminants" annotation (Dialog(enable = false));

  .Modelica.Units.SI.MoleFraction y_cont[N, nCont] "Contamint species molar fractions in each node";

protected 

  final parameter Integer iCont[nCont] = Medium_an.substanceIndexVector(contaminantsSpecies) "Index of considered contaminants in the anodic medium" annotation (Dialog(enable = false));
 
initial equation  
  
  for i in 1:nCont loop
   
    assert(iCont[i] > 0, "Membrane: the selected medium for the anode side does not contain the required substances", level = AssertionLevel.error);
   
  end for;

equation
  
  for i in 1:nCont loop

    y_cont[:, i] = y_an[:, iCont[i]];

  end for;

end BaseContaminantsLoss;