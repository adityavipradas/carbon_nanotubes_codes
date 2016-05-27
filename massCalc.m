function[mass] = massCalc(CNT_num, CNT_density, CNT_diameter, CNT_thick, CNT_length)
%CNT_density gm/cc to kg/m3
CNT_density = CNT_density * 1000;
%CNT_diameter nm to m
CNT_diameter = CNT_diameter * (10^-9);
%CNT_thick pm to m
CNT_thick = CNT_thick * (10^-12);
%CNT_length um to m
CNT_length = CNT_length * (10^-6);
volume = (pi/4) * (CNT_diameter^2 - (CNT_diameter - (2 * CNT_thick))^2) * CNT_length;
mass = (volume * CNT_density) * 10^6
sprintf('Mass of %d CNTs is %f mg',CNT_num, mass)
