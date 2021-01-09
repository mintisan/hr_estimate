function [currentState, logP] = hmmviterbi(tr,e)

numStates = size(tr,1);
checkE = size(e,1);
numEmissions = size(e,2);

% work in log space to avoid numerical issues
L = numEmissions;
currentState = zeros(1,L);
if L == 0
    return
end
w = warning('off'); % get log of zero warnings
logTR = log(tr);
logE =  (e);
%keyboard;
warning(w);

% allocate space
pTR = zeros(numStates,L);

% !!!!!!!!!!!!!!!!!!!!!!!!!!!!
% assumption is that model is in the most likely state at step 0
%v(1,1) = 0;
%vOld = v;

v = repmat(-inf, numStates,1);
vOld = diag(logTR).*logE(:,1);


% loop through the model
for count = 1:L
    for state = 1:numStates
        % for each state we calculate
        % v(state) = e(state,seq(count))* max_k(vOld(:)*tr(k,state));
        bestVal = -inf;
        bestPTR = 0;
        % use a loop to avoid lots of calls to max
        for inner = 1:numStates 
            val = vOld(inner) + logTR(inner,state);
            if val > bestVal
                bestVal = val;
                bestPTR = inner;
            end
        end
        % save the best transition information for later backtracking
        pTR(state,count) = bestPTR;
        % update v
        v(state) = logE(state,count) + bestVal;
    end
    vOld = v;
end

% decide which of the final states is post probable
[logP, finalState] = max(v);

% Now back trace through the model
currentState(L) = finalState;
for count = L-1:-1:1
    currentState(count) = pTR(currentState(count+1),count+1);
    if currentState(count) == 0
        error('stats:hmmviterbi:ZeroTransitionProbability',...
       ['A zero transition probability was encountered from state %d.\n',...
       'Please provide more data or PseudoTransition information.'],...
       currentState(count+1));
    end
end


