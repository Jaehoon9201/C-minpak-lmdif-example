x  = [-1.7237128,1.8712276,-0.96608055,...
      -0.28394297,1.3416969,1.3757038,...
      -1.3703436,0.042581975,-0.14970151,...
       0.82065094];

y     = [0.19000429,6.5807428,1.4582725,...
      2.7270851,5.5969253,5.6249280,...
      0.787615,3.2599759,2.9771762,...
      4.5936475];


error_factor = 1;
a0 = 3.18970;
a1 = 1 ;
a2 = 0.149020 ; 

opt = a0 + a1 .* x + a2.* x.* x ;

plot(y, 'k'); hold on;
plot(opt, 'b');
plot(error_factor.*(y-opt), 'r--'); grid on;
legend('real vals', 'opt vals', 'residuals');
