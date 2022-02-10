function M = energy_two(Nos,Nsq,ak2,L)

M = eye(Nos+Nsq,Nos+Nsq);
Cos = two_two(Nos,L);
Dos = deven_two(Nos,L);
Wos = (Dos'*Cos*Dos+ak2*Cos);
Wsq = two_two(Nsq,L);
% Wos = (Dos'*Cos*Dos+ak2*Cos)/(ak2);
% Wsq = two_two(Nsq)/ak2;

[u,s,v]=svd(Wos); s=sqrt(diag(s));
Mos=diag(s)*u';

[u,s,v]=svd(Wsq); s=sqrt(diag(s));
Msq = diag(s)*u';

M = [Mos zeros(Nos,Nsq); zeros(Nsq,Nos) Msq];

end