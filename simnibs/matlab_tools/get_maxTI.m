function absEmax = get_maxTI(E1,E2)
% absEmax = get_maxTI(E1,E2)
% 
% calculates the maximal modulation amplitude of the TI envelope using 
% the equation given in Grossman et al, Cell 169, 1029–1041.e6, 2017
%
% E1 and E2 have to be Nx3 vectors, or element_data structures containing
% tridata and tetdata that are Nx3 vectors
%
% A. Thielscher, 2020

if ~isstruct(E1)
    absEmax=get_maxTI_nx3(E1,E2);
else
    absEmax = struct('name','maxTIamplitude','tridata',[],'tetdata',[]);
    absEmax.tridata = get_maxTI_nx3(E1.tridata,E2.tridata);
    absEmax.tetdata = get_maxTI_nx3(E1.tetdata,E2.tetdata);
end


function absEmax = get_maxTI_nx3(E1,E2)
    % calculate amplitude of low-frequency envelope for two Nx3 vectors
    
    if (size(E1,2)~=3)||(size(E2,2)~=3); error('both vectors have to be Nx3'); end
    if size(E1,2)~=size(E2,2); error('both vectors have to have the same size'); end

    % ensure E1>E2
    idx=sum(E2.^2,2)>sum(E1.^2,2);
    E1hlp=E1;
    E1(idx,:)=E2(idx,:);
    E2(idx,:)=E1hlp(idx,:);

    % ensure alpha < pi/2
    idx=dot(E1,E2,2)<0;
    E2(idx,:)=-E2(idx,:);

    % get |Emax|
    absE1=sqrt(sum(E1.^2,2));
    absE2=sqrt(sum(E2.^2,2));
    cosalpha=dot(E1,E2,2)./(absE1.*absE2);

    absEmax=2*sqrt(sum(cross(E2,E1-E2,2).^2,2))./sqrt(sum((E1-E2).^2,2));
    idx=absE2<=absE1.*cosalpha;
    absEmax(idx)=2*absE2(idx);




