

%% K-Middle Mean Evaluation
function mean=KMM(k,Vector)
Lenght_of_Vector = numel(Vector);
Vector=sort(Vector);
if mod(Lenght_of_Vector,2) == 1
    half_lengh = 0.5*(Lenght_of_Vector+1);
    factor = 1/(2*k - 1);
    K_middle_vector=Vector((half_lengh-k+1):(half_lengh+k-1));
    sum_of_element=sum(K_middle_vector);
else
    half_lengh=0.5*Lenght_of_Vector;
    factor=1/(2*k);
    K_middle_vector=Vector((half_lengh-k+1):(half_lengh+k));
    sum_of_element=sum(K_middle_vector);
end
mean=factor*sum_of_element;
end
















