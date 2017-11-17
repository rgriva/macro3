% Função auxiliar INV_DIST_NUM para calcular distribuição invariante da matriz
% de Markov M, cujas linhas somam 1, pelo metodo iterativo.
% M eh a matriz em que estamos interessados e epsilon eh o parâmetro de 
% convergência. 

function inv_dist = inv_dist_num(M, epsilon)
inv_dist0 = ones(length(M),1);
inv_dist1 = M'*inv_dist0;
diff = norm(inv_dist1 - inv_dist0);

while diff > epsilon
    inv_dist0 = inv_dist1;
    inv_dist1 = M'*inv_dist0;
    diff = norm(inv_dist1 - inv_dist0);
end

inv_dist = inv_dist1/sum(inv_dist1);