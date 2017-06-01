% This file calls get_subbands.m to get the parameters dependent on N, and
% then prints the parameters in subbands.txt in a format that can be used
% in gnrfet.lib

% Change the numbers below to specify the range of N
N_start = 6;
N_end = 50;


% DO NOT MODIFY ANYTHING BELOW
n_sbbd=3;

fid = fopen('subbands.txt','w');

first_line=1;
for N0=N_start:N_end
    N=N0;
    get_subbands_3NN();
    if (first_line==1)
        fprintf(fid, '.if (N == %d)\n', N);
        first_line=0;
    else
        fprintf(fid, '.elseif (N == %d)\n', N);
    end
    for s=1:n_sbbd
        fprintf(fid, '    .param sbbd%da = ''%d''\n', s, sbbd(s));
    end
    for s=1:n_sbbd
        fprintf(fid, '    .param mass%da = ''%d''\n', s, mass(s));
    end
%     fprintf(fid, '    .param sbbd2a = ''%d''\n', sbbd2);
%     fprintf(fid, '    .param sbbd3a = ''%d''\n', sbbd3);
%     fprintf(fid, '    .param mass1a = ''%d''\n', Mass1);
%     fprintf(fid, '    .param mass2a = ''%d''\n', Mass2);
%     fprintf(fid, '    .param mass3a = ''%d''\n', Mass3);
    N=N0-1;
    get_subbands_3NN();
    for s=1:n_sbbd
        fprintf(fid, '    .param sbbd%db = ''%d''\n', s, sbbd(s));
    end
    for s=1:n_sbbd
        fprintf(fid, '    .param mass%db = ''%d''\n', s, mass(s));
    end
%     fprintf(fid, '    .param sbbd1b = ''%d''\n', sbbd1);
%     fprintf(fid, '    .param sbbd2b = ''%d''\n', sbbd2);
%     fprintf(fid, '    .param sbbd3b = ''%d''\n', sbbd3);
%     fprintf(fid, '    .param mass1b = ''%d''\n', Mass1);
%     fprintf(fid, '    .param mass2b = ''%d''\n', Mass2);
%     fprintf(fid, '    .param mass3b = ''%d''\n', Mass3);
end
fprintf(fid, '.endif\n');
fclose(fid);