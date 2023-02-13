%% Metabolite:metabolite square matrix of shared reactions 

function [metmet, metmet2] = metmetmat(rxnList)


for m1 = 1:length(rxnList)
    for m2 = 1:length(rxnList)
        mat1 = rxnList{m1}{1,2};
        mat2 = rxnList{m2}{1,2};

        if m1 == m2
            continue 
        else 
            if any(ismember(mat1, mat2))
                metmet{m1, m2} = intersect(mat1,mat2); 
                metmet2(m1,m2) = 1; 
            end 
        end 
    end
end



end 