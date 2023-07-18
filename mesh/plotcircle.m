function plotcircle(varargin)

hold on

if nargin==3
    x0 = varargin{1};  y0 = varargin{2}; r = varargin{3};
end
if nargin==2
    p0 = varargin{1};  r = varargin{2};
    x0 = p0(1); y0 = p0(2);
end

t = linspace(0,2*pi,720);
x = x0 + r*cos(t);
y = y0 + r*sin(t);

plot(x,y,'r','linewidth',1);

hold off

