
clear all;
% syms pi;
syms tau;
syms tauc;

s = 3; %packet size (slots)

i = 0; %counter stage 
j = 5; %contention window size (slots)
k = 10; %ATIM window size (slots)
nn = 4; %nodes number, including myself

% k = w+s+1; %for computing purpose 
%(i,j,k) represents one state in Markov Chain



A=zeros(j,k*nn); %column, row
A = sym(A);

B=zeros(j,(s+1)*nn); %column, row
B = sym(B);



% A(:,:,2)=0;
% A(:,:,3)=0;

% for n=1:j*k
%     
%     A(n) = ps+n;
%     
% end
% 
% 
% disp(A(5,14));

for q=j:-1:1
    At=zeros(j,k*nn); %column, row
    At = sym(At);
    At(q,1) = 1;
    
    Bt=zeros(j,(s+1)*nn); %column, row
    Bt = sym(Bt);
    
    for m=q:-1:1
        for n=1:k*nn
%             disp(m);
%             disp(n);
            if At(m,n)~=0
                if m>1 
%                     disp('m>1');
                    no = nn - mod(n,nn);
                    if(no==nn)
                        no=0;
                    end
%                     disp(no);
                    
                    if n+1*nn<=k*nn
                        At(m-1,n+1*nn) = At(m,n)*tauc^no+At(m-1,n+1*nn);
%                         if(no==0)
%                             fprintf('A(%d,%d) is same as A(%d,%d)\n', m,n,m-1,n+1*nn);
%                         end
                    end          
                    if n+s*nn<=k*nn && no~=0
                        At(m-1,n+s*nn+1) = At(m,n)*(nchoosek(no,1)*tau*tauc^(no-1))+At(m-1,n+s*nn);
                        At(m-1,n+s*nn) = At(m,n)*(1-tauc^no-nchoosek(no,1)*tau*tauc^(no-1))+At(m-1,n+s*nn);
                    end
                    
                    if n+1*nn>k*nn
                        Bt(m-1,n+1*nn-k*nn) = At(m,n)*tauc^no+Bt(m-1,n+1*nn-k*nn);
%                         if(no==0)
%                             fprintf('A(%d,%d) is same as B(%d,%d)\n', m,n,m-1,n+1*nn-k*nn);
%                             disp(At(m,n));
%                             disp(Bt(m-1,n+1*nn-k*nn));
%                         end
%                         disp(Bt(m-1,n+1-k));
                    end

                    if n+s*nn>k*nn && no~=0
                        Bt(m-1,n+s*nn-k*nn+1) = At(m,n)*(nchoosek(no,1)*tau*tauc^(no-1))+Bt(m-1,n+s*nn-k*nn);
                        Bt(m-1,n+s*nn-k*nn) = At(m,n)*(1-tauc^no-nchoosek(no,1)*tau*tauc^(no-1))+Bt(m-1,n+s*nn-k*nn);
%                         disp(Bt(m-1,n+s-k));
                    end
                 end
            end
        end
    end
%     At.'
%     Bt.'
    A=A+At;
    B=B+Bt;
end



pr=zeros(j,1); %column, row
pr = sym(pr);

for q=1:j
        
        for n=1:(s+1)*nn
            pr(q,1) = B(q,n) + pr(q,1);
        end
    
end

A = A.';%simple transpose, not conjugate transpose
disp('A:');
disp(A);

B = B.';%simple transpose, not conjugate transpose
disp('B:');
disp(B);

disp('pr:');
pr = pr.';
disp(pr);

tau = 0.2;
tauc = 1 - tau;
% ps = (1-tau)^8;

for q=1:j
        
    disp(eval(pr(q)));
    
end

% 
% if A(1)=='2'
%     disp('ffffff');
% else
%     disp('dddddd');
% end


% syms a;
% 
% b = a*ps;
% c = b/a;





% 
% disp(c);