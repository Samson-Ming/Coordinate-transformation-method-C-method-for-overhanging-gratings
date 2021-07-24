switch grating
    case 0
        Polyline1
    case 1
        Polyline2    
    case 2
        Full_Triangle         
    case 3
        Full_Trapezoid
    case 4
        Harmonics1
    case 5
        Harmonics2 
    case 6
        Sampling      
    case 7
        Power_sine
    case 8
        %Better put the definition file in the grating profiles filefolder. 
        eval(input_grating_file(1:end-2))   
end