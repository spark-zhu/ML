function [aAndB aInds bInds] = d_uIntersect(A,B)

%
%
%
% Don Vaughn

[orderedAAndB orderedAInds orderedBInds] = intersect(A,B);
[aInds sortToUnsortedMap]= sort(orderedAInds);
bInds = orderedBInds(sortToUnsortedMap);
aAndB = A(aInds);