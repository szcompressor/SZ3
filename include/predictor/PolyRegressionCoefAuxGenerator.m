% Matlab code to compute the inverse matrix of poly coefficient aux matrix
N=10;
BLOCK_SIZE_MIN=3;
BLOCK_SIZE_MAX=10;
% generate_coef_inverse(N,6,6,6);

fileID = fopen('PolyRegressionCoefAux.bin','w');

for i=BLOCK_SIZE_MIN:BLOCK_SIZE_MAX
    for j=BLOCK_SIZE_MIN:BLOCK_SIZE_MAX
        for k=BLOCK_SIZE_MIN:BLOCK_SIZE_MAX
            Ai=generate_coef_inverse(N,i,j,k);
            fwrite(fileID,Ai,'single');
        end
    end
end
fclose(fileID);

function Ai = generate_coef_inverse(N, n1, n2,n3)
    A=zeros(N);
    for i = 0:n1-1
        for j=0:n2-1
            for k=0:n3-1
                single_row=[1, i, j, k, i * i, i * j, i * k, j * j, j * k, k * k];
                for p1 = 1:N
                    for p2=1:N
                        A(p1,p2)= A(p1, p2)+single_row(p1)*single_row(p2);
                    end
                end
            end
        end
    end
    Ai=inv(A);
%     disp(A);
%     disp(Ai);
end