function g = gamma_phi(phases)
%%%
%%% g = gamma_phi(phases)
%%%
%%% generates the gamma matrix/operator used to convert 
%%% a vector (like power injection) to a phase x phase matrix.
%%% phases is a vector of the present phases. for example [1,3] would
%%% mean that only phases 1 and 3 are present (i.e. phase 2 is absent)

a = exp(-1i*2*pi/3);
gamma = [1  , a^2, a  ;
         a  , 1  , a^2;
         a^2, a  , 1];
g = gamma(phases,phases);
