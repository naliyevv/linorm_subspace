function [X,Y,Z] = plotcontours(A,B,xmin,xmax,ymin,ymax)


hx = 1;
hy = 0.05;
n = size(A);

j = 1;
for x = xmin:hx:xmax
	k = 1;
	for y = ymin:hy:ymax
		X(k,j) = x;
		Y(k,j) = y;
		Z(k,j) = min(svd([A-(x+i*y)*eye(n)  B]));
		k = k+1;
	end

	j = j+1;
end

return;