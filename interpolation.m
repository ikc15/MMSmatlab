function out= interpolation(in,n)
% in = input vector 
% n = generate n sample points. 
% https://uk.mathworks.com/help/matlab/ref/interp1.html
x = 1:length(in);
xq = linspace(1,length(in),n);
out = interp1(x,in,xq,'spline');
% plot(x,in,'o',xq,out,':.');
% legend('original', 'interpolated');
% title('Spline interpolation')
end