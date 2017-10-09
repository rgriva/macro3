% Lista 2 - Macroeconomia II (2017)
% Aluno: Raul Guarini Riva
% Questão 3 - item 3
% Função auxiliar para calcular distribuição invariante da matriz
% estocástica M, cujas linhas somam 1, numericamente.
% M eh a matriz em que estamos interessados e epsilon eh o parâmetro de convergência. 

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