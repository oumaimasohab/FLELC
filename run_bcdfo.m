clear all

for nprob = [101]

   %  Initialize test problem from file obj.m

   if ( (nprob >= 1 & nprob <= 71) || (nprob >=101 & nprob <= 171) )

      [ x0, xl, xu, nprob ] = bcdfo_init( nprob );

   else

      disp(['Problem ', num2str(nprob), ' does not exist in bcdfo_init().']);
      return

   end

   %  Possibility to modify default values of parameters

   tolmg      = 1e-5;      % stopping criterion for gradient norm
   maxeval    = 1000;      % stopping criterion for nbr. of f evaluations
   whichmodel = 0;         % choice of local interpolation models
   hardcons   = 0;         % respect bounds as hard constraints
   noisy      = 0;         % function is supposed to be noisy

   [ x, fx, gx, msg, nit, nf ] = bcdfo( @(x) objf(x,nprob), x0, ...
      'lbounds', xl, 'ubounds', xu, ...
      'epsilon', tolmg, ...
      'maxeval', maxeval, 'maxit', maxeval, ...
      'whichmodel', whichmodel, ...
      'hardcons', hardcons, ...
      'noisy', noisy);

nprob
end
