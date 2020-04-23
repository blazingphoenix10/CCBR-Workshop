 function [total] = makeSparsity(amps)                     
   
    total = 1-((nansum(amps).^2)./(nansum(amps.^2)))/length(amps);

 end