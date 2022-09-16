function [ trel ] = cctrellis(constraintLength,genMatrix)

K = size(genMatrix,1); 
N = size(genMatrix,2); %(n,k,m) code
trel.numInputSymbols = 2^K;
trel.numOutputSymbols = 2^N;
numReg = sum(constraintLength) - K;%m
trel.numStates = 2^numReg; % m^2 states
%A trellis is represented by a structure with the following fields:

%numInputSymbols – Number of input symbols

%numOutputSymbols – Number of output symbols

%numStates – Number of states

%nextStates – Next state matrix


trel.nextStates = zeros(trel.numStates,trel.numInputSymbols);
outputsN = zeros(trel.numStates,trel.numInputSymbols,N); %outputs – Output matrix
nStates = zeros(trel.numStates,trel.numInputSymbols);
currentStates = 0:2^numReg-1; 
multiFactor = 0;
for k = 1:K
    a = 2^(constraintLength(k)-1);
    currentStatesSub = mod(currentStates,a);
    currentStatesShift = floor(currentStatesSub/2);
    currentStates = floor(currentStates/a); %return the nearest integer less than 
    % or equal to that currentStates/a.

    index_1 = (mod(floor((0:2^K-1)/2^(K-k)),2) == 1); %next states
    index_0 = (mod(floor((0:2^K-1)/2^(K-k)),2) == 0); 
    nStates(:,index_0) = repmat(currentStatesShift',1,2^(K-1));
    nStates(:,index_1) = repmat(currentStatesShift',1,2^(K-1)) + a/2;
    trel.nextStates = trel.nextStates + ...
        2^(multiFactor)*nStates; %update the states
    multiFactor = multiFactor + constraintLength(k) -1;
    
    currentStatesSubBit = dec2bin(currentStatesSub) - 48; %creat the ouput
    for n = 1:N
        genBit = dec2bin(base2dec(num2str(genMatrix(k,n)),8),constraintLength(k)) - 48;
        rOutput = currentStatesSubBit.*repmat(genBit(2:end),2^numReg,1);
        rOutput = mod(sum(rOutput,2),2);
        rOutput = repmat(rOutput,1,2^K);
        rOutput(:,index_1) =  rOutput(:,index_1) + genBit(1);
        outputsN(:,:,n)  = mod(outputsN(:,:,n) + rOutput,2);
    end
end
trel.outputs = zeros(trel.numStates,trel.numInputSymbols);
for n = 1:N
    trel.outputs = trel.outputs + outputsN(:,:,n)*2^(N-n); % update the output trellis
end
end

