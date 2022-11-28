
function out=delay_finder_starter_2(A,fs,binning)

del=zeros(size(A,1),size(A,1));
start=find(A(:,1)==1);
for kk1=1:size(A,1)
    if isempty(find(A(kk1,2:end),1,'first'))==0 && ismember(kk1,start)==0
        del(start,kk1)=find(A(kk1,2:end),1,'first');
    end
end   
out=del.*((1/fs)*binning);

end