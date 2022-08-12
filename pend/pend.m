function dxdt = pend(t,x,uex)

  global T_samp
  
  % ZOH the control vector
  try
      u = uex(fix(t/T_samp)+1);
  catch
      u = uex(end);
  end
  
  g = 9.81;
  l = 0.5;
  m = 0.15;
  v = 0.1;
  
  dx1 = x(2);
  dx2 = (g/l)*sin(x(1)) - (v/(m*l^2))*x(2) + (1/(m*l^2))*u;

  dxdt = [dx1; dx2];
  
end