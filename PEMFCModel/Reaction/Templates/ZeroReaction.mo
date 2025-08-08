within PEMFCModel.Reaction.Templates;
model ZeroReaction
  
  extends .PEMFCModel.Reaction.Templates.DynamicReaction(final S_reac = zeros(0, Medium.nS), final kEq_reac = ones(N, n_reac));

end ZeroReaction;