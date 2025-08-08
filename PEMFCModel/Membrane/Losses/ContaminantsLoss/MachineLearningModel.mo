within PEMFCModel.Membrane.Losses.ContaminantsLoss;

model MachineLearningModel "Voltage loss due to contaminants by means of a machine learning model"
 
// Model to be implemented!
  
  extends PEMFCModel.Membrane.Losses.Templates.BaseContaminantsLoss;

  parameter Real k = 10e-2;
    
equation

  for i in 1:N loop
  
    E_loss_cell[i] = k*(y_cont[i,:]*y_cont[i,:]);
    
 end for;
    
 annotation(Icon(graphics={Text(
          extent={{-64,30},{68,-28}},
          textColor={28,108,200},
          textString="E_cont")}));

end MachineLearningModel;