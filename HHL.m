clear all
close all
clc

e0 = [1;0]; e1 = [0;1];
e00 = [1;0;0;0]; e01 = [0;1;0;0]; e10 = [0;0;1;0]; e11 = [0;0;0;1];
H = [1,1;1,-1]/sqrt(2);

% randomly generate matrix A in HHL
[U D] = svd(rand(4));
u1 = U(:,1); u2 = U(:,2); u3 = U(:,3); u4 = U(:,4);
l = ceil(rand(1,4)*2);
l = ((l - 1) .*2) + 1;
l1 = l(1); l2 = l(2); l3 = l(3); l4 = l(4);
A = l1*u1*u1'+l2*u2*u2'+l3*u3*u3'+l4*u4*u4';
t0 = 2*pi;
T = 4;

% randomly generate vector b in HHL
[Ub D] = svd(rand(4));
b = Ub(:,1);
%Ux: |0> -> |x>
x = inv(A)*b; x = x/norm(x);
Ux = [x,null(x*x')];
%Uf: p,q
Uf = kron(e00*e00',expm(1i*A*0*t0/T)) + kron(e01*e01',expm(1i*A*1*t0/T)) + kron(e10*e10',expm(1i*A*2*t0/T)) + kron(e11*e11',expm(1i*A*3*t0/T));
QFT = [1,1,1,1;1,1i,-1,-1i;1,-1,1,-1;1,-1i,-1,1i]/2;
QFTinv = QFT';
C = 0.7;
Uc1 = kron(e0,e00); 
Uc2 = kron(sqrt(1-C^2/1^2)*e0+C/1*e1,e01);
Uc3 = kron(sqrt(1-C^2/2^2)*e0+C/2*e1,e10);
Uc4 = kron(sqrt(1-C^2/3^2)*e0+C/3*e1,e11);
Ucg = [Uc1,Uc2,Uc3,Uc4];
Ucgn = null(Ucg*Ucg');
Uc = [Ucg,Ucgn];    %r, p
M0 = e1*e1'; M1 = e0*e0'; Ir = eye(2); Iq = eye(4); Ip = eye(4); Ia = eye(2);

% the projection operators for the assertions
%implement R: First Measure p1, p2: result: 0,0
%Introduce auxiliary qubit a = |0>, then apply Ux' on q, Us on a,r,q
%Measure a: result 0; Recover: apply Us' on a,r,q, apply Ux on q.
Us = eye(16); Us(6:11,6:11) = fliplr(eye(6));
R = kron(kron(e1*e1',e00*e00'),x*x')+kron(kron(e0*e0',e00*e00'),Iq);
%implement P: Measure r, p1, p2: result: 0,0,0
P = kron(kron(e0*e0',e00*e00'),Iq);
%implement Q: Apply Ux' on q, then Measure r, p1, p2, q1, q2:
%result:1,0,0,0,0, then Apply Ux on q.
Q = kron(kron(e1*e1',e00*e00'),x*x');
%implement S:
S = kron(Ir, kron(e01 * e01' + e11 * e11', Iq));

% begin execute the HHL program
disp("begin HHL")
%
state = kron(kron(e0,e00),e00);

% first time enter loop
state = kron(M1, kron(Ip, Iq)) * state;

% implement assertion P
disp("assertion P")
disp("1st time entering the loop")
state = kron(kron(e0*e0',e00*e00'),Iq) * state; % apply the measurement in assertion P
observe_Probability_Distribution = state'.'.*state;
disp(observe_Probability_Distribution(1:4))

state = kron(kron(Ir,Ip),Ub)*state;
state = kron(kron(Ir,kron(H,H)),Iq)*state;
state = kron(Ir,Uf)*state;
state = kron(kron(Ir,QFTinv),Iq)*state;

% implement assertion S
disp("assertion S")
disp("before assertion S")
observe_Probability_Distribution = state'.'.*state;
disp(observe_Probability_Distribution(5:8))
disp(observe_Probability_Distribution(13:16))
disp("after assertion S")
state = kron(Ir, kron(e01 * e01' + e11 * e11', Iq)) * state; % apply assertion S
observe_Probability_Distribution = state'.'.*state;
disp(observe_Probability_Distribution(5:8))
disp(observe_Probability_Distribution(13:16))

state = kron(Uc,Iq)*state;
state = kron(kron(Ir,QFT),Iq)*state;
state = kron(Ir,Uf')*state;
state = kron(kron(Ir,kron(H,H)),Iq)*state;

% implement assertion R
disp("assertion R")
disp("before entering assertion R")
observe_Probability_Distribution = state'.'.*state;
disp(observe_Probability_Distribution(1:4))
disp(observe_Probability_Distribution(17:20))
disp("after the introduced unitary transformations Ux and Ur")
state = kron(Ir, kron(Ip, Ux')) * state; % apply Ux
state = kron(e0, state); % introduce anxiliary qubit a
% Ur is omited since it will not affect the state vector
observe_Probability_Distribution = state'.'.*state;
disp(observe_Probability_Distribution(1:4))
disp(observe_Probability_Distribution(17:20))
% Ur' is omited since it will not affect the state vector
state = state(1:32); %  remove anxiliary qubit a
state = kron(Ir, kron(Ip, Ux)) * state; % apply Ux'

% store the current state
temp_state = state;

% apply the measurement and exit the while loop
state = kron(M0,kron(Ip,Iq))*state;  
state = state/sqrt((state'*state));  % normalize the state vector


% implement assertion Q
disp("assertion Q")
disp("before entering assertion Q")
observe_Probability_Distribution = state'.'.*state;
disp(observe_Probability_Distribution(17:20))
disp("after the introduced unitary transformation Ux")
state = kron(Ir, kron(Ip, Ux')) * state;
observe_Probability_Distribution = state'.'.*state;
disp(observe_Probability_Distribution(17:20))
state = kron(Ir, kron(Ip, Ux)) * state;

% restore state before the last measurement
state = temp_state;

% apply the measurement and re-enter the while loop
state = kron(M1,kron(Ip,Iq))*state;  
state = state/sqrt((state'*state));  % normalize the state

% implement assertion P
disp("assertion P")
disp("2nd time entering the loop")
state = kron(kron(e0*e0',e00*e00'),Iq) * state; % apply the measurement in assertion P
observe_Probability_Distribution = state'.'.*state;
disp(observe_Probability_Distribution(1:4))

%state(17:20),
%state(17:20)./x
%x
