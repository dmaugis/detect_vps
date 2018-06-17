function ortho_vps = orthogonal_triplet(my_vps, NFAs, ORTHOGONALITY_THRESHOLD);
% returns most significant orthogonal triplet
N = size(my_vps,2);


nfa_scores_triplets = [];
ortho_scores_triplets =[];
triplets = [];


nfa_scores_pairs = [];
ortho_scores_pairs =[];
pairs = [];

for i=1:N
    vpi = my_vps(:,i);
    nfai = NFAs(i);
    
    for j=1:i-1
        if i == j
            continue
        end
        
        vpj = my_vps(:,j);
        nfaj = NFAs(j);
        scoreij = abs(dot(vpi,vpj));
        
        nfa_scores_pairs = [nfa_scores_pairs nfai+nfaj];
        ortho_scores_pairs = [ortho_scores_pairs scoreij];
        pairs = [pairs; i j];
        
        
        for k=1:N
            if k==j || j==i
                continue
            end
            
            
            vpk = my_vps(:,k);
            nfak = NFAs(k);
            scorejk = abs(dot(vpj,vpk));
            scoreik = abs(dot(vpi,vpk));
            
            % get orthogonality score
            
            ortho_score = max([ scoreij,  scorejk,  scoreik]);
            nfa_score = nfai+nfaj+nfak;
            
            nfa_scores_triplets = [nfa_scores_triplets nfa_score];
            ortho_scores_triplets = [ortho_scores_triplets ortho_score];
            triplets = [triplets; i j k];
            
        end
    end
end

ORTHO_SCORE_THRESHOLD_TRIPLETS = ORTHOGONALITY_THRESHOLD;
ORTHO_SCORE_THRESHOLD_pairS = ORTHOGONALITY_THRESHOLD;

z3 = find(ortho_scores_triplets > ORTHO_SCORE_THRESHOLD_TRIPLETS);
nfa_scores_triplets(z3)=[];
ortho_scores_triplets(z3)=[];
triplets(z3,:)=[];


z2 = find(ortho_scores_pairs > ORTHO_SCORE_THRESHOLD_pairS);

nfa_scores_pairs_orig = nfa_scores_pairs;
ortho_scores_pairs_orig = ortho_scores_pairs;
pairs_orig = pairs;

nfa_scores_pairs(z2)=[];
ortho_scores_pairs(z2)=[];
pairs(z2,:)=[];


if isempty(triplets)
    % no orthogonal triplets, return pair
    
     if isempty(pairs)
%         [V, I] = sort(NFAs, 'descend');
%         ortho_vps = my_vps(:,I(1));
        
        
        % if no triplets or pairs, get the most orthogonal pair
        [V I] = sort(ortho_scores_pairs_orig);
        pair = pairs_orig(I(1),:);
        ortho_vps = [my_vps(:,pair(1)) my_vps(:,pair(2))];
        
    else
        % by nfa
        [V, I] = sort(nfa_scores_pairs, 'descend');
        % by orthogonality
        %[V, I] = sort(ortho_scores_pairs, 'ascend');
        
        pair = pairs(I(1),:);
        ortho_vps = [my_vps(:,pair(1)) my_vps(:,pair(2))];
    end
else
    % there are triplets
    
    % by nfa
    [V, I] = sort(nfa_scores_triplets, 'descend');
    % by orthogonnality
    %[V, I] = sort(ortho_scores_triplets, 'ascend');
    
    
    triplet = triplets(I(1),:);
    ortho_vps = [my_vps(:,triplet(1)) my_vps(:,triplet(2)) my_vps(:,triplet(3))];
end


