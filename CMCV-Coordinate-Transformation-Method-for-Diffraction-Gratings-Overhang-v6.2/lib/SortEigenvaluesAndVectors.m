function [eig_p,vect_p,eig_m,vect_m]=SortEigenvaluesAndVectors(eig,vect,beta_m_Pos, beta_m_Neg,accuracyImag)
%% sortation of eigenvalues and corresp eigenvectros
eig_real_p=[];
vec_real_p=[];
eig_real_m=[];
vec_real_m=[];

if ~isempty(beta_m_Pos)
    idxReal=find(abs(imag(eig)) < accuracyImag);
    eig_rot=eig(idxReal);
    idxRealAndPositive=[];
    idxRealAndNegative=[];
    for iii=1:numel(beta_m_Pos)
        [~,index]=min(abs(eig_rot-beta_m_Pos(iii)));
        idxRealAndPositive=[idxRealAndPositive,idxReal(index)];
    end

     for iii=1:numel(beta_m_Neg)
        [~,index]=min(abs(eig_rot-beta_m_Neg(iii)));
        idxRealAndNegative=[idxRealAndNegative,idxReal(index)];
 end
%in fact not neccessary
eig_real_p = eig(idxRealAndPositive).';
vec_real_p = vect(:,idxRealAndPositive);
eig_real_m = eig(idxRealAndNegative).';
vec_real_m = vect(:,idxRealAndNegative);
end

%length(eig_real_p)
%length(eig_real_m)

idxImagAndPositive = abs(imag(eig)) >= accuracyImag & imag(eig) > 0;
idxImagAndNegative = abs(imag(eig)) >= accuracyImag & imag(eig) <= 0;


eig_comp_p = eig(idxImagAndPositive).';
%length(eig_comp_p)
vec_comp_p = vect(:,idxImagAndPositive);
eig_comp_m = eig(idxImagAndNegative).';
%length(eig_comp_m)
vec_comp_m = vect(:,idxImagAndNegative);

%% sortation in descending order
if ~isempty(eig_real_p)
    [~,ind]=sort(real(eig_real_p),'descend');
    eig_real_p=eig_real_p(ind);
    vec_real_p=vec_real_p(:,ind);
else
    eig_real_p =ones(1,0);
    vec_real_p =ones(size(vect,1),0);
end

if ~isempty(eig_comp_p)
    [~,ind]=sort(imag(eig_comp_p),'ascend');
    eig_comp_p=eig_comp_p(ind);
    vec_comp_p=vec_comp_p(:,ind);
end

if ~isempty(eig_real_m)
    [~,ind]=sort(real(eig_real_m),'ascend');
    eig_real_m=eig_real_m(ind);
    vec_real_m=vec_real_m(:,ind);
else
    eig_real_m =ones(1,0);
    vec_real_m =ones(size(vect,1),0);
end

if ~isempty(eig_comp_m)
    [~,ind]=sort(imag(eig_comp_m),'descend');
    eig_comp_m=eig_comp_m(ind);
    vec_comp_m=vec_comp_m(:,ind);
end

eig_p = [eig_real_p, eig_comp_p];
vect_p = [vec_real_p, vec_comp_p];
eig_m = [eig_real_m, eig_comp_m];
vect_m = [vec_real_m, vec_comp_m];

end