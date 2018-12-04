function [L,M,N]=HybridMinVarAnalysis(fgm, B, t1, t2)
% Hybrid Minimum Variance Analysis

% find index of init and end time 
ix1 = Index(fgm{1,1},t1); % use the Index.m function
ix2 = Index(fgm{1,1},t2);

% check interval selected correctly:
%disp('Intervals between:')
%disp(datetime(fgm{1,1}(ix1),'ConvertFrom','datenum','Format','dd-MMM-yyyy HH:mm:ss.SSS'));
%disp(datetime(fgm{1,1}(ix2),'ConvertFrom','datenum','Format','dd-MMM-yyyy HH:mm:ss.SSS'));

% perform MVA on B-field in the interval
[out,l,v]=irf_minvar(B(ix1:ix2, 1:3));
% l3 is smallest eigenvalue- corresponds to N
% l1 is largest eigenvalue- corresponds to L
% v = row vectors?? manual says v is column vectors but writes:
% v(1,:)=first vector
% perform hybrid minimum variance
B1 = mean(B(ix1-1:ix1,1:3));
B2 = mean(B(ix2:ix2+1,1:3));
N = cross(B1,B2)/norm(cross(B1,B2));
M = cross(N,v(1,:));
L = cross(M,N);

end