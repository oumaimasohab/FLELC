%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Testing Script - Validation CUTEst test problems
%
% This script calls the CUTEst/MATLAB interface to check if the problems as
% listed have correct dimensions.
%
% In the same directory, it is assumed that there exists a file ListPbms 
% containing the data related to every problem in the collection under the 
% form:
%
% NameProblem	Dimension	BestValue
%
% An output file is produced under the following form
%
% NameProblem	Dimensions	InitialValue	InitialGradientNorm	...
% InitialMinimumEigenvalue
%
% Implementation: C. W. Royer 
% Last updated: 2020.07.27
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
%
% Path to the CUTEst/MATLAB interface (Example on Clement's laptop)
addpath('/home/ouma/cutest/src/matlab')

% Getting the problem list
fid = fopen('ListProblems','r');
PROBS = textscan(fid,'%s %d %f');
fclose(fid);
% Opening the output file
fid2 = fopen('ResProblems','w');
%
npbs=length(PROBS{1});
DIMS = PROBS{2};
BVALS = PROBS{3};
%
%npbs=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Testing the extraction of the problem data from CUTEst
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
for numpb = 1:npbs
%
%	Decoding the SIF file of the problem
   	st=(['!cutest2matlab ' PROBS{1}{numpb}]);
   	disp(st); eval(st); pause(1);
%
%	Gathering the problem data
   	prob = cutest_setup();
   	x0=prob.x;
	fbest = BVALS(numpb);
	n= prob.n;
   	name=prob.name;
%
	if (DIMS(numpb)~=n)
		error('CUTEST:pb %s has wrong dimension',name);
	end
%
	fprintf('%s & %d\n',name,n);
	f0 = cutest_obj(x0);
	g0 = cutest_grad(x0);
	if norm(x0)~=0
		l0 = (x0'*cutest_hprod(x0,x0))/(x0'*x0);
	else
		l0 = 0;
	end
	fprintf(fid2,'%s & %d &%1.3e & %1.3e & %1.3e\n',name,n,f0,norm(g0),l0);
%	pause;
%
end
%
fclose(fid2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%