function Final_distance = Electrode_Unit_distance(Location, Boundary)
% Boundary(1) is the left most side x corrdinate. 
% Boundary(2) is the right most side x corrdinate.
% Boundary(3) is the lower most side y corrdinate. 
% Boundary(4) is the upper most side y corrdinate.
% Location(1) is the x corrdinate of the unit 
% Location(2) is the y corrdinate of the unit



distance(1) = Boundary(1)-Location(1);
distance(2) = -Boundary(2)+Location(1);
distance(3) = Boundary(3)-Location(2);
distance(4) = -Boundary(4)+Location(2);

Final_distance = max(distance);




