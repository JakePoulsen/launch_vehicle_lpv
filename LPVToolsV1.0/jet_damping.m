function J = jet_damping(Theta, Mass_flow, Mass_dir, Center_mass, Inertia, Mass)

%k is the axis radius of gyration [3x1]

%Mass is the mass flow rate of mass out the nozzles in kg/second. [1x1]

%Center_mass is the distance from the center of gravity of the rocket to
% - the exit of the nozzles [1x1]

%Theta is the angular velocity of the given axis. [3x1]

%Mass_dir is the direction of the mass outflow(same direction as thrust [3x1]


k = (I/Mass)^(0.5);
%moment of inertia
%area cross section

Mass_dir_norm = abs((Mass_dir)/(norm(Mass_dir)));
Mass_x = Mass_flow*(Mass_dir_norm(2)+Mass_dir_norm(3));
Mass_y = Mass_flow*(Mass_dir_norm(1)+Mass_dir_norm(3));
Mass_z = Mass_flow*(Mass_dir_norm(1)+Mass_dir_norm(2));
J_x = -Mass_x*Theta(1)*(Center_mass^2 - k(1)^2);
J_y = -Mass_y*Theta(2)*(Center_mass^2 - k(2)^2);
J_z = -Mass_z*Theta(3)*(Center_mass^2 - k(3)^2);

J = [J_x; J_y; J_z]; %The return values
        return;
end