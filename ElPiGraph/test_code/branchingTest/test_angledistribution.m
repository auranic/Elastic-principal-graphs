npoints = 1000;
dimension = 10;

data = zeros(npoints,dimension);

for i=1:npoints
  r = (rand()-0.5)*2;
  %r = rand();
  var = exp(-(r)*(r)*3)/5;
  %var = 0;
  noise = 0.01;
  round = exp(-r*r*10);
  div = 10;
  if r<0
    data(i,1) = -r+round/div+var*(rand()-0.5);
    data(i,2) = round/div+var*(rand()-0.5)+rand()*noise;
    data(i,3) = var*(rand()-0.5)+rand()*noise;
  else
    data(i,2) = r+round/div+var*(rand()-0.5);
    data(i,1) = round/div+var*(rand()-0.5)+rand()*noise;
    data(i,3) = var*(rand()-0.5)+rand()*noise;
  end
  for k=3:dimension
      data(i,k) = var*(rand()-0.5)/10+rand()*noise;
  end
end

for i=npoints+1:npoints+npoints/3
    r = rand();
    var = exp(-(r)*(r)*3)/10;
    data(i,1) = r+0.2+var*(rand()-0.5)+rand()*noise;
    data(i,2) = -r*r+var*(rand()-0.5)+rand()*noise;
    data(i,3) = -r*r+var*(rand()-0.5)+rand()*noise;
    for k=4:dimension
      data(i,k) = rand()*noise;
    end
end
    