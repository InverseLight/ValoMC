function distance = distanceFromLine(p1,p2,p3)
% calculates the distance from point p3 from a line that goes trough points p1 and p2
  x0 = p3(1);
  y0 = p3(2);
  x1 = p1(1);
  y1 = p1(2);
  x2 = p2(1);
  y2 = p2(2);
  distance = abs(((y2-y1)*x0 - (x2-x1)*y0 +x2*y1- y2*x1))/norm(p2-p1);
end
