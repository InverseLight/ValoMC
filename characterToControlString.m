function c = characterToControlString(str, mode)
   c = {};

   if(mode == 1)
       % lightsources
       for k = 1:length(str)
           if(str(k) == 'a') c(k) = {'none'}; end;
           if(str(k) == 'l') c(k) = {'direct'}; end;
           if(str(k) == 'c') 
               c(k) = {'cosinic'}; 
           end;
           if(str(k) == 'g') c(k) = {'gaussian'}; end;
           if(str(k) == 'i') c(k) = {'isotropic'}; end;
           if(str(k) == 'p') c(k) = {'pencil'}; end;                        
       end
       if(mode == 2)
           % light source directions
           if(str(k) == 'a') c(k) = {'absolute'}; end;
           if(str(k) == 'r') c(k) = {'relative'}; end;
       end
   end

end

