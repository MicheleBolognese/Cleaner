within PEMFCModel.Membrane.MassTransport.WaterContent;
function WaterContentGandomiFunction 
"It returns the water content as function of water activity according to Gandomi's correlation"
   
  extends .Modelica.Icons.Function;

  input Real a "Water activity";
    
  output Real lambda "Water content";
    
  algorithm
    
    lambda := Modelon.Math.Smoothing.spliceFunction(16.8, Modelon.Math.Smoothing.spliceFunction(14 + 1.4*(a - 1), 0.043 + 17.18*a - 39.85*a^2 + 36.0*a^3, a-1, 0.01), a-3, 0.01);
    
end WaterContentGandomiFunction;
