function d = EuclideanDistance(a,b)
% DISTANCE - computes Euclidean distance matrix
%
% E = EuclideanDistance(A,B)
%
%    A - (MxD) matrix 
%    B - (NxD) matrix
%
% Returns:
%    E - (MxN) Euclidean distances between vectors in A and B
%
%
% Description : 
%    This fully vectorized (VERY FAST!) m-file computes the 
%    Euclidean distance between two vectors by:
%
%                 ||A-B|| = sqrt ( ||A||^2 + ||B||^2 - 2*A.B )
%
% Example : 
%    A = rand(100,400); B = rand(200,400);
%    d = EuclideanDistance(A,B);

aa=sum(a.^2,2); bb=sum(b.^2,2); ab=a*b';
d = sqrt(-bsxfun(@minus,bsxfun(@minus,2*ab, aa), bb')); 
%d = sqrt(abs(repmat(aa,[1 size(bb,1)]) + repmat(bb',[size(aa,1) 1]) - 2*ab));