function [bcattrib] = PLUG_bcfunction_con(flagptr,a);
global numcase bcflagc centelem bedge coord;

switch numcase
    case 231
        bcattrib=bcflagc(flagptr,2);
    case 232
        bcattrib=bcflagc(flagptr,2);
    case 233
        bcattrib=bcflagc(flagptr,2);
    case 234
        bcattrib=bcflagc(flagptr,2);
    case 235
       bcattrib=bcflagc(flagptr,2); 
    case 236
       bcattrib=bcflagc(flagptr,2); 
    case 237
       bcattrib=bcflagc(flagptr,2);
    case 238
       bcattrib=bcflagc(flagptr,2);
    case 239
       bcattrib=bcflagc(flagptr,2); 
    case 241
       bcattrib=bcflagc(flagptr,2);  
    case 242
       bcattrib=bcflagc(flagptr,2);
    case 243
       bcattrib=bcflagc(flagptr,2); 
    case 244
       bcattrib=bcflagc(flagptr,2); 
    case  245
      bcattrib=bcflagc(flagptr,2);
    case  247
      bcattrib=bcflagc(flagptr,2);
     case  249
      bcattrib=bcflagc(flagptr,2);
    case 248
        
        %The initial condition is a integral mean of function.
        
        bcattrib = sin(0.25*pi*(flagptr+2*a));
      
      
    case 246
      bcattrib=bcflagc(flagptr,2);  
end
end