%--------------------------------------------------------------------
% PURPOSE:
%
% This routine assembles the local matrices corresponding
% to the neighbours of a given element.
%
% Author:
% Paola Antonietti
%--------------------------------------------------------------------

function [M]= assemble_neigh_A(M,row,neight,M1,nln,n_edge,nA) %nP

for iedg=1:n_edge
    nP=0;
    if (neight(iedg) > 0  && neight(iedg) >= nP + 1 && neight(iedg) <= nA)
        neight_A = neight(iedg) - nP;
        j=(neight_A-1)*nln*ones(nln,1) + [1:nln]';
        M(row,j)=M(row,j)+M1(:,:,iedg);
    end
end