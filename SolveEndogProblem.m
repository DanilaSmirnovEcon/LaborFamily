% This function solves the single household problem
function solutionSP = SolveEndogProblem(solutionSfguess, solutionSmguess, loadsol, dispiter, param, gridSf, gridSm, gridP)
always_marry = param.always_marry;
tremble = param.tremble;
param = ensureCBlocks(param);

function param = ensureCBlocks(param)
    if ~isfield(param,'CJO'), param.CJO = param.CJU; end
    if ~isfield(param,'COJ'), param.COJ = param.CUJ; end
    if ~isfield(param,'COU'), param.COU = param.CUU; end
    if ~isfield(param,'CUO'), param.CUO = param.CUU; end
    if ~isfield(param,'COO'), param.COO = param.CUU; end
end

%interpolation indecies from single big grid to pair smaller grid

fVai_m = zeros(gridSf.a.n,2);
fVai_m_w = zeros(gridSf.a.n,2);
for ai_ = 1:gridSf.a.n
    [ai_m_diff, ai_m_diffi] = sort(abs(gridSf.a.a(ai_) - gridP.af.a));
    fVai_m(ai_,:) = ai_m_diffi(1:2);
    den = ai_m_diff(1) + ai_m_diff(2);
    if den == 0
        fVai_m_w(ai_,:) = [1, 0];  % exact hit on first node
    else
        fVai_m_w(ai_,:) = [ai_m_diff(2)/den, ai_m_diff(1)/den];
    end
end

fVwi_m = zeros(gridSf.w.n,2);
fVwi_m_w = zeros(gridSf.w.n,2);
for wi_ = 1:gridSf.w.n
    [wi_m_diff, wi_m_diffi] = sort(abs(gridSf.w.w(wi_) - gridP.wf.w));
    fVwi_m(wi_,:) = wi_m_diffi(1:2);
    den = wi_m_diff(1) + wi_m_diff(2);
    if den == 0
        fVwi_m_w(wi_,:) = [1, 0];
    else
        fVwi_m_w(wi_,:) = [wi_m_diff(2)/den, wi_m_diff(1)/den];
    end
end

mVai_m = zeros(gridSm.a.n,2);
mVai_m_w = zeros(gridSm.a.n,2);
for ai_ = 1:gridSm.a.n
    [ai_m_diff, ai_m_diffi] = sort(abs(gridSm.a.a(ai_) - gridP.am.a));
    mVai_m(ai_,:) = ai_m_diffi(1:2);
    den = ai_m_diff(1) + ai_m_diff(2);
    if den == 0
        mVai_m_w(ai_,:) = [1, 0];
    else
        mVai_m_w(ai_,:) = [ai_m_diff(2)/den, ai_m_diff(1)/den];
    end
end

mVwi_m = zeros(gridSm.w.n,2);
mVwi_m_w = zeros(gridSm.w.n,2);
for wi_ = 1:gridSm.w.n
    [wi_m_diff, wi_m_diffi] = sort(abs(gridSm.w.w(wi_) - gridP.wm.w));
    mVwi_m(wi_,:) = wi_m_diffi(1:2);
    den = wi_m_diff(1) + wi_m_diff(2);
    if den == 0
        mVwi_m_w(wi_,:) = [1, 0];
    else
        mVwi_m_w(wi_,:) = [wi_m_diff(2)/den, wi_m_diff(1)/den];
    end
end


uJl_f = instUtilityJl(gridSf,param.pf);
uUl_f = instUtilityUl(gridSf,param.pf);
uLl_f = [uJl_f(:); uUl_f(:)];
uJc_f = instUtilityJc(gridSf,param.pf);
uUc_f = instUtilityUc(gridSf,param.pf);
uLc_f = [uJc_f(:); uUc_f(:)];

uOl_f = instUtilityOl(gridSf,param.pf);
uOc_f = instUtilityOc(gridSf,param.pf);


uJl_m = instUtilityJl(gridSm,param.pm);
uUl_m = instUtilityUl(gridSm,param.pm);
uLl_m = [uJl_m(:); uUl_m(:)];
uJc_m = instUtilityJc(gridSm,param.pm);
uUc_m = instUtilityUc(gridSm,param.pm);
uLc_m = [uJc_m(:); uUc_m(:)];

uOl_m = instUtilityOl(gridSm,param.pm);
uOc_m = instUtilityOc(gridSm,param.pm);



uJJl = instUtilityJJl();
uJUl = instUtilityJUl();
uJOl = instUtilityJOl();
uUJl = instUtilityUJl();
uOJl = instUtilityOJl();
uUUl = instUtilityUUl();
uOUl = instUtilityOUl();
uUOl = instUtilityUOl();
uOOl = instUtilityOOl();

uJJc = instUtilityJJc();
uJUc = instUtilityJUc();
uJOc = instUtilityJOc();
uUJc = instUtilityUJc();
uOJc = instUtilityOJc();
uUUc = instUtilityUUc();
uUOc = instUtilityUOc();
uOUc = instUtilityOUc();
uOOc = instUtilityOOc();
uLLl = [reshape(uJJl, param.NP,1);...
    reshape(uJUl, param.NP/gridP.wm.n,1);...
    reshape(uUJl, param.NP/gridP.wf.n,1);...
    reshape(uUUl, gridP.af.n*gridP.am.n,1)];
uLLc = [reshape(uJJc, param.NP,1);...
    reshape(uJUc, param.NP/gridP.wm.n,1);...
    reshape(uUJc, param.NP/gridP.wf.n,1);...
    reshape(uUUc, gridP.af.n*gridP.am.n,1)];
 

uOLl = [reshape(uOJl, param.NP/gridP.wf.n,1);...
    reshape(uOUl, gridP.af.n*gridP.am.n,1)];
uOLc = [reshape(uOJc, param.NP/gridP.wf.n,1);...
    reshape(uOUc, gridP.af.n*gridP.am.n,1)];

uLOl = [reshape(uJOl, param.NP/gridP.wm.n,1);...
    reshape(uUOl, gridP.af.n*gridP.am.n,1)];
uLOc = [reshape(uJOc, param.NP/gridP.wm.n,1);...
    reshape(uUOc, gridP.af.n*gridP.am.n,1)];
 


try
    if loadsol
        load('solutionSP','solutionSP');
    else
        error('no load');
    end
catch
    % Initialize the values guess
    solutionSP.JflValue = solutionSfguess.JlValue;
    solutionSP.JmlValue = solutionSmguess.JlValue;
    solutionSP.UflValue = solutionSfguess.UlValue;
    solutionSP.UmlValue = solutionSmguess.UlValue;
    solutionSP.OflValue = solutionSfguess.OlValue;
    solutionSP.OmlValue = solutionSmguess.OlValue;

    solutionSP.JJlValue = uJJl;
    solutionSP.JUlValue = uJUl;
    solutionSP.UJlValue = uUJl;
    solutionSP.UUlValue = uUUl;
    solutionSP.efl=solutionSfguess.el./2;
    solutionSP.eml=solutionSmguess.el./2;
    solutionSP.hfl=solutionSfguess.hl./2;
    solutionSP.hml=solutionSmguess.hl./2;
    solutionSP.ofl=solutionSfguess.ol./2;
    solutionSP.oml=solutionSmguess.ol./2;

    solutionSP.eel=abs(uJJl);
    solutionSP.ehl=abs(uJUl);
    solutionSP.hel=abs(uUJl);
    solutionSP.hhl=abs(uUUl);
    solutionSP.eol=abs(uJOl);
    solutionSP.oel=abs(uOJl);
    solutionSP.hol=abs(uUOl);
    solutionSP.ohl=abs(uOUl);
    solutionSP.ool=abs(uOOl);

    solutionSP.JfcValue = solutionSfguess.JcValue;
    solutionSP.JmcValue = solutionSmguess.JcValue;
    solutionSP.UfcValue = solutionSfguess.UcValue;
    solutionSP.UmcValue = solutionSmguess.UcValue;
    solutionSP.OfcValue = solutionSfguess.OcValue;
    solutionSP.OmcValue = solutionSmguess.OcValue;

    solutionSP.JJcValue = uJJc;
    solutionSP.JUcValue = uJUc;
    solutionSP.UJcValue = uUJc;
    solutionSP.UUcValue = uUUc;

    solutionSP.efc=solutionSfguess.ec./2;
    solutionSP.emc=solutionSmguess.ec./2;
    solutionSP.hfc=solutionSfguess.hc./2;
    solutionSP.hmc=solutionSmguess.hc./2;
    solutionSP.ofc=solutionSfguess.oc./2;
    solutionSP.omc=solutionSmguess.oc./2;

    solutionSP.eec=abs(uJJc);
    solutionSP.ehc=abs(uJUc);
    solutionSP.hec=abs(uUJc);
    solutionSP.hhc=abs(uUUc);
    solutionSP.eoc=abs(uJOc);
    solutionSP.oec=abs(uOJc);
    solutionSP.hoc=abs(uUOc);
    solutionSP.ohc=abs(uOUc);
    solutionSP.ooc=abs(uOOc);

    solutionSP.OJlValue = uOJl;
    solutionSP.OUlValue = uOUl;
    solutionSP.OJcValue = uOJc;
    solutionSP.OUcValue = uOUc;
    
    solutionSP.JOlValue = uJOl;
    solutionSP.UOlValue = uUOl;
    solutionSP.JOcValue = uJOc;
    solutionSP.UOcValue = uUOc;

    solutionSP.OOlValue = uOOl;
    solutionSP.OOcValue = uOOc;


    % Define the Act vectors
    solutionSP.ActLl = ones(param.pf.NSall+param.pm.NSall,1);
    solutionSP.ActLLl_LOl = ones(param.NP+param.NP/gridP.wm.n+param.NP/gridP.wf.n+gridP.af.n*gridP.am.n,1);
    solutionSP.ActLLl_OLl = ones(param.NP+param.NP/gridP.wf.n+param.NP/gridP.wm.n+gridP.af.n*gridP.am.n,1);
    solutionSP.ActLc = solutionSP.ActLl;
    solutionSP.ActLLc_LOc = solutionSP.ActLLl_LOl;
    solutionSP.ActLLc_OLc = solutionSP.ActLLl_OLl;
    
    solutionSP.ActOLl_OOl = ones(param.NP/gridP.wf.n+gridP.af.n*gridP.am.n,1);
    solutionSP.ActOLc_OOc = solutionSP.ActOLl_OOl;
    
    solutionSP.ActLOl_OOl = ones(param.NP/gridP.wm.n+gridP.af.n*gridP.am.n,1);
    solutionSP.ActLOc_OOc = solutionSP.ActLOl_OOl;
    
    % Define the U_O vectors
    solutionSP.Lfl_Ofl = 1*ones(gridSf.a.n,1);
    solutionSP.OLl_LLl = 1*ones(param.NP/gridP.wf.n+gridP.af.n*gridP.am.n,1);
    solutionSP.Lfc_Ofc = solutionSP.Lfl_Ofl ;
    solutionSP.OLc_LLc = solutionSP.OLl_LLl;
    
    solutionSP.Lml_Oml = 1*ones(gridSm.a.n,1);
    solutionSP.LOl_LLl = 1*ones(param.NP/gridP.wm.n+gridP.af.n*gridP.am.n,1);
    solutionSP.Lmc_Omc = solutionSP.Lml_Oml;
    solutionSP.LOc_LLc = solutionSP.LOl_LLl;
    
    solutionSP.OO_LO = 1*ones(2*gridP.af.n*gridP.am.n,1);
    solutionSP.OO_OL = solutionSP.OO_LO;

    solutionSP.constr = -Inf;

end

% 
uOOl = [reshape(uOOl, gridP.af.n*gridP.am.n,1)];
uOOc = [reshape(uOOc, gridP.af.n*gridP.am.n,1)];

% Solve for the values
difff=1; diffS=1;

% Marriage-market pool: the stationary distribution of singles in the
% no-marriage problem, solved with the same parameters and grids as the
% full problem. It is held fixed (not updated with the equilibrium) so the
% pool of potential partners does not become an extra fixed point.
solutionSf_s = SolveSingleProblem(param, param.pf, gridSf);
solutionSm_s = SolveSingleProblem(param, param.pm, gridSm);
solutionSP.solutionSf_s = solutionSf_s;
solutionSP.solutionSm_s = solutionSm_s;


nJf = gridSf.a.n*gridSf.w.n;
nJm = gridSm.a.n*gridSm.w.n;
nOf = gridSf.a.n; nUf=nOf;
nOm = gridSm.a.n; nUm=nOm;
nJJ = gridP.af.n*gridP.wf.n*gridP.am.n*gridP.wm.n; %JJ size
nJO = gridP.af.n*gridP.wf.n*gridP.am.n; nJU=nJO;%JU size
nOJ = gridP.af.n*gridP.am.n*gridP.wm.n; nUJ=nOJ;%UJ size
nOO = gridP.af.n*gridP.am.n; nUU=nOO;nOU=nOO;nUO=nOO;%UU size

nALLL = 2*(nJf+nUf+nJm+nUm+nJJ+nJU+nUJ+nUU); %2* because of no kid (l) and kid (c)
nAOOL = 2*(nOf+nOJ+nOU);
nAOLO = 2*(nOm+nJO+nUO);
nAOO  = 2*(nOO);

for ji =1:param.maxiter 
    
    [solutionSP, diffS, diffSOm, diffSOf, diffPJJ, diffPJO, diffPOJ, diffPOO] = updateValue(solutionSP,solutionSf_s,solutionSm_s,param, gridP, gridSf, gridSm);
    
    
    if dispiter
        disp("iter "+num2str(ji)+...
            "; diffS "+num2str(diffS)+...
            "; diffSOm "+num2str(diffSOm)+...
            "; diffSOf "+num2str(diffSOf)+...
            "; diffPJJ "+num2str(diffPJJ)+...
            "; diffPJO "+num2str(diffPJO)+...
            "; diffPOJ "+num2str(diffPOJ)+...
            "; diffPOO "+num2str(diffPOO)+...
            "; difff "+num2str(difff))
    end
    if (max([diffS,diffPJJ])<param.critC) && (ji>1)
        if dispiter
        disp("Endog problem converged in "+num2str(ji)+...
            " iterations; diff distribution "+num2str(difff)+...
            " ; diff value "+num2str([diffS, diffSOm, diffSOf, diffPJJ, diffPJO, diffPOJ, diffPOO]))
        end
        break
    end
end
if diffS>param.critC
    disp(diffS);
    disp("WARNING!!! Agent value function not found");   
    save('solutionSP_partial','solutionSP','param','gridSf','gridSm','gridP');
    %error("err")
end

% updating 
function [solutionSP, diffL, diffSOm, diffSOf, diffPLL, diffPLO, diffPOL, diffPOO] = updateValue(solutionSP,solutionSf_s,solutionSm_s,param, gridP, gridSf, gridSm)
    
    %%===================================================================
    % Solve problem for the labor force participants
    %%===================================================================
       
    %%= Genetrate transition matrixes    
    % tic;
    % single part for female
    [AfJl_Jl, AfJl_Ul, AfJl_Ol, AfJl_JJl, AfJl_JUl, AfJl_JOl, AfJl_JJc, AfJl_JUc, AfJl_JOc] = generateAJl_f(solutionSP, solutionSm_s,param.pf, gridSf);
    [AfUl_Jl, AfUl_Ul, AfUl_Ol, AfUl_UJl, AfUl_UUl, AfUl_UOl, AfUl_UJc, AfUl_UUc, AfUl_UOc] = generateAUl_f(solutionSP, solutionSm_s,param.pf, gridSf);
    [AfJc_Jc, AfJc_Uc, AfJc_Oc, AfJc_JJc, AfJc_JUc, AfJc_JOc] = generateAJc_f(solutionSP, solutionSm_s,param.pf, gridSf);
    [AfUc_Jc, AfUc_Uc, AfUc_Oc, AfUc_UJc, AfUc_UUc, AfUc_UOc] = generateAUc_f(solutionSP, solutionSm_s,param.pf, gridSf);        
    AfLl = [AfJl_Jl, AfJl_Ul;...
            AfUl_Jl, AfUl_Ul];
    AfLc = [AfJc_Jc, AfJc_Uc;...
            AfUc_Jc, AfUc_Uc];    
    AfLl = AfLl - param.pf.pcs.*speye(param.pf.NSall,param.pf.NSall);
    % % female single (no kids): Jf->Of and Jf->JO ; Uf->Of and Uf->UO
    % Df_LLl = blkdiag( ...
    %     spdiags(sum(AfJl_Ol,2) + sum(AfJl_JOl,2), 0, nJf, nJf), ...
    %     spdiags(sum(AfUl_Ol,2) + sum(AfUl_UOl,2), 0, nUf, nUf) );
    % AfLl = AfLl - Df_LLl;
    % 
    % % female single (kids): Jf^c->Of^c and Jf^c->JO^c ; Uf^c->Of^c and Uf^c->UO^c
    % Df_LLc = blkdiag( ...
    %     spdiags(sum(AfJc_Oc,2) + sum(AfJc_JOc,2), 0, nJf, nJf), ...
    %     spdiags(sum(AfUc_Oc,2) + sum(AfUc_UOc,2), 0, nUf, nUf) );
    % AfLc = AfLc - Df_LLc;

    % single part for male
    [AmJl_Jl, AmJl_Ul, AmJl_Ol, AmJl_JJl, AmJl_UJl, AmJl_OJl, AmJl_JJc, AmJl_UJc, AmJl_OJc] = generateAJl_m(solutionSP, solutionSf_s,param.pm, gridSm);
    [AmUl_Jl, AmUl_Ul, AmUl_Ol, AmUl_JUl, AmUl_UUl, AmUl_OUl, AmUl_JUc, AmUl_UUc, AmUl_OUc] = generateAUl_m(solutionSP, solutionSf_s,param.pm, gridSm);
    [AmJc_Jc, AmJc_Uc, AmJc_Oc, AmJc_JJc, AmJc_UJc, AmJc_OJc] = generateAJc_m(solutionSP, solutionSf_s,param.pm, gridSm);
    [AmUc_Jc, AmUc_Uc, AmUc_Oc, AmUc_JUc, AmUc_UUc, AmUc_OUc] = generateAUc_m(solutionSP, solutionSf_s,param.pm, gridSm);        
    AmLl = [AmJl_Jl, AmJl_Ul;...
            AmUl_Jl, AmUl_Ul];
    AmLc = [AmJc_Jc, AmJc_Uc;...
            AmUc_Jc, AmUc_Uc];    
    AmLl = AmLl - param.pm.pcs.*speye(param.pm.NSall,param.pm.NSall); 
    % % male single (no kids): Jm->Om and Jm->OJ ; Um->Om and Um->OU
    % Dm_LLl = blkdiag( ...
    %     spdiags(sum(AmJl_Ol,2) + sum(AmJl_OJl,2), 0, nJm, nJm), ...
    %     spdiags(sum(AmUl_Ol,2) + sum(AmUl_OUl,2), 0, nUm, nUm) );
    % AmLl = AmLl - Dm_LLl;
    % 
    % % male single (kids)
    % Dm_LLc = blkdiag( ...
    %     spdiags(sum(AmJc_Oc,2) + sum(AmJc_OJc,2), 0, nJm, nJm), ...
    %     spdiags(sum(AmUc_Oc,2) + sum(AmUc_OUc,2), 0, nUm, nUm) );
    % AmLc = AmLc - Dm_LLc;
        
    % double part

    [AJJl_JJl, AJJl_JUl, AJJl_UJl, AJJl_JOl, AJJl_OJl] = generateAJJl(solutionSP, param, gridP);
    [AJUl_JUl, AJUl_JJl, AJUl_UUl, AJUl_OUl, AJUl_JOl] = generateAJUl(solutionSP, param, gridP);
    [AUJl_UJl, AUJl_JJl, AUJl_UUl, AUJl_OJl, AUJl_UOl] = generateAUJl(solutionSP, param, gridP); 
    [AUUl_UUl, AUUl_UJl, AUUl_JUl, AUUl_UOl, AUUl_OUl] = generateAUUl(solutionSP, param, gridP);

    [AJJc_JJc, AJJc_JUc, AJJc_UJc, AJJc_JOc, AJJc_OJc] = generateAJJc(solutionSP, param, gridP);
    [AJUc_JUc, AJUc_JJc, AJUc_UUc, AJUc_OUc, AJUc_JOc] = generateAJUc(solutionSP, param, gridP);
    [AUJc_UJc, AUJc_JJc, AUJc_UUc, AUJc_OJc, AUJc_UOc] = generateAUJc(solutionSP, param, gridP);    
    [AUUc_UUc, AUUc_UJc, AUUc_JUc, AUUc_UOc, AUUc_OUc] = generateAUUc(solutionSP, param, gridP);
            
    % Transitions for double
    AJJ_UU=sparse(zeros(nJJ,nUU));
    AJU_UJ=sparse(zeros(nJU,nUJ));
    AUJ_JU=sparse(zeros(nUJ,nJU));
    AUU_JJ=sparse(zeros(nUU,nJJ));
    ALLl = [AJJl_JJl, AJJl_JUl, AJJl_UJl, AJJ_UU;...
            AJUl_JJl, AJUl_JUl, AJU_UJ,   AJUl_UUl;...
            AUJl_JJl, AUJ_JU,   AUJl_UJl, AUJl_UUl;...
            AUU_JJ,   AUUl_JUl, AUUl_UJl, AUUl_UUl];
    ALLl = ALLl - param.pcp.*speye(param.NPall,param.NPall); 
    % D_LL_toO_l = blkdiag( ...
    % spdiags(sum(AJJl_OJl,2) + sum(AJJl_JOl,2), 0, nJJ, nJJ), ...
    % spdiags(sum(AJUl_OUl,2) + sum(AJUl_JOl,2), 0, nJU, nJU), ...
    % spdiags(sum(AUJl_OJl,2) + sum(AUJl_UOl,2), 0, nUJ, nUJ), ...
    % spdiags(sum(AUUl_OUl,2) + sum(AUUl_UOl,2), 0, nUU, nUU) );
    % ALLl = ALLl - D_LL_toO_l; % last two check

    ALLc = [AJJc_JJc, AJJc_JUc, AJJc_UJc, AJJ_UU;...
            AJUc_JJc, AJUc_JUc, AJU_UJ,   AJUc_UUc;...
            AUJc_JJc, AUJ_JU,   AUJc_UJc, AUJc_UUc;...
            AUU_JJ,   AUUc_JUc, AUUc_UJc, AUUc_UUc];
    % D_LL_toO_c = blkdiag( ...
    % spdiags(sum(AJJc_OJc,2) + sum(AJJc_JOc,2), 0, nJJ, nJJ), ...
    % spdiags(sum(AJUc_OUc,2) + sum(AJUc_JOc,2), 0, nJU, nJU), ...
    % spdiags(sum(AUJc_OJc,2) + sum(AUJc_UOc,2), 0, nUJ, nUJ), ...
    % spdiags(sum(AUUc_OUc,2) + sum(AUUc_UOc,2), 0, nUU, nUU) );
    % ALLc = ALLc - D_LL_toO_c; % last two check

    % Transitions from single female to double
    AfJ_UJ = zeros(nJf,nUJ);
    AfJ_UU = zeros(nJf,nUU);
    AfU_JJ = zeros(nUf,nJJ);
    AfU_JU = zeros(nUf,nJU);
    AfLl_LLl = [AfJl_JJl, AfJl_JUl, AfJ_UJ,   AfJ_UU;...
                AfU_JJ,   AfU_JU,   AfUl_UJl, AfUl_UUl];
    AfLl_LLl = sparse(AfLl_LLl);
    AfLl_LLc = [AfJl_JJc, AfJl_JUc, AfJ_UJ,   AfJ_UU;...                        
                AfU_JJ,   AfU_JU,   AfUl_UJc, AfUl_UUc];
    AfLl_LLc = sparse(AfLl_LLc);
    AfLc_LLc = [AfJc_JJc, AfJc_JUc, AfJ_UJ,   AfJ_UU;...                        
                AfU_JJ,   AfU_JU,   AfUc_UJc, AfUc_UUc];
    AfLc_LLc = sparse(AfLc_LLc);

    % Transitions from single male to double
    AmJ_JU = zeros(nJm,nJU);
    AmJ_UU = zeros(nJm,nUU);
    AmU_JJ = zeros(nUm,nJJ);
    AmU_UJ = zeros(nUm,nUJ);
    AmLl_LLl = [AmJl_JJl, AmJ_JU,   AmJl_UJl, AmJ_UU;...
                AmU_JJ,   AmUl_JUl, AmU_UJ,   AmUl_UUl];
    AmLl_LLl = sparse(AmLl_LLl);    
    AmLl_LLc = [AmJl_JJc, AmJ_JU,   AmJl_UJc, AmJ_UU;...
                AmU_JJ,   AmUl_JUc, AmU_UJ,   AmUl_UUc];
    AmLl_LLc = sparse(AmLl_LLc);
    AmLc_LLc = [AmJc_JJc, AmJ_JU,   AmJc_UJc, AmJ_UU;...
                AmU_JJ,   AmUc_JUc, AmU_UJ,   AmUc_UUc];
    AmLc_LLc = sparse(AmLc_LLc);

    % Build aggregate matrix
    ALLLl = sparse([AfLl, zeros(size(AfLl,1),size(AmLl,1)), AfLl_LLl;...
                zeros(size(AmLl,1),size(AfLl,1)), AmLl, AmLl_LLl;...
                    zeros(size(ALLl,1), size(AfLl,1)),zeros(size(ALLl,1), size(AmLl,1)), ALLl]);
    ALLLc = sparse([AfLc, zeros(size(AfLc,1),size(AmLc,1)), AfLc_LLc;...
                zeros(size(AmLc,1),size(AfLc,1)), AmLc, AmLc_LLc;...
                    zeros(size(ALLc,1), size(AfLc,1)),zeros(size(ALLc,1), size(AmLc,1)), ALLc]);
    ALLLl_LLLc = [ ...
              param.pf.pcs.*speye(param.pf.NSall), zeros(param.pf.NSall, param.pm.NSall), AfLl_LLc; 
              zeros(param.pm.NSall, param.pf.NSall), param.pm.pcs.*speye(param.pm.NSall), AmLl_LLc; 
              sparse(param.NPall, param.pf.NSall),  sparse(param.NPall, param.pm.NSall),  param.pcp.*speye(param.NPall) ...
            ];
    ALLLc_LLLl = sparse(param.NPall+param.pf.NSall+param.pm.NSall,param.NPall+param.pf.NSall+param.pm.NSall);
    ALLL = [ALLLl     , ALLLl_LLLc;...
            ALLLc_LLLl, ALLLc];
    solutionSP.ALLL = ALLL;
    
          
    u_LLL = [uLl_f;uLl_m;uLLl;uLc_f;uLc_m;uLLc];


    %%===================================================================
    % Solve problem for the labor force Female nonparticipant
    %%===================================================================
    [AfOl_Ol, AfOl_Jl, AfOl_OJl, AfOl_OUl, AfOl_OJc, AfOl_OUc] = generateAOl_f(solutionSP, solutionSm_s, param.pf, gridSf);
    [AfOc_Oc, AfOc_Jc, AfOc_OJc, AfOc_OUc] = generateAOc_f(solutionSP, solutionSm_s, param.pf, gridSf);
    AfOl = AfOl_Ol;
    AfOc = AfOc_Oc;    
    AfOl = AfOl - param.pf.pcs.*speye(gridSf.a.n,gridSf.a.n);
    % ---- FEMALE O(single) ----
    % AfOl = AfOl_Ol ...
    %      - param.pf.pcs .* speye(nOf) ...
    %      - spdiags(sum(AfOl_Jl,2), 0, nOf, nOf);      % Ol -> Jl outflow
    % 
    % AfOc = AfOc_Oc ...
    %      - spdiags(sum(AfOc_Jc,2), 0, nOf, nOf);      % Oc -> Jc outflow


    [AOJl_OJl, AOJl_OUl, AOJl_OOl, AOJl_JJl] = generateAOJl(solutionSP, param, gridP);
    [AOUl_OUl, AOUl_OJl, AOUl_OOl, AOUl_JUl] = generateAOUl(solutionSP, param, gridP);
    [AOJc_OJc, AOJc_OUc, AOJc_OOc, AOJc_JJc] = generateAOJc(solutionSP, param, gridP);
    [AOUc_OUc, AOUc_OJc, AOUc_OOc, AOUc_JUc] = generateAOUc(solutionSP, param, gridP);

    AOLl = [AOJl_OJl, AOJl_OUl;...
         AOUl_OJl, AOUl_OUl];
    AOLl = AOLl - param.pcp.*speye(param.NP/gridP.wf.n+gridP.af.n*gridP.am.n,param.NP/gridP.wf.n+gridP.af.n*gridP.am.n);
    AOLc = [AOJc_OJc, AOJc_OUc;...
         AOUc_OJc, AOUc_OUc];
    % ---- COUPLE O (OJ, OU; no kids) ----
    % AOLl = [AOJl_OJl, AOJl_OUl;
    %         AOUl_OJl, AOUl_OUl] ...
    %      - param.pcp .* speye(nOJ+nOU);
    % AOLl = AOLl ...
    %      - blkdiag( spdiags(sum(AOJl_JJl,2),0,nOJ,nOJ), ... % OJ_l -> JJ_l
    %                 spdiags(sum(AOUl_JUl,2),0,nOU,nOU) );    % OU_l -> JU_l
    % AOLc = [AOJc_OJc, AOJc_OUc;
    %         AOUc_OJc, AOUc_OUc];
    % AOLc = AOLc ...
    %      - blkdiag( spdiags(sum(AOJc_JJc,2),0,nOJ,nOJ), ... % OJ_c -> JJ_c
    %                 spdiags(sum(AOUc_JUc,2),0,nOU,nOU) );    % OU_c -> JU_c
    
    AfOl_OLl = [AfOl_OJl, AfOl_OUl];
    AfOl_OLl = sparse(AfOl_OLl);
    AfOl_OLc = [AfOl_OJc, AfOl_OUc];
    AfOl_OLc = sparse(AfOl_OLc);
    AfOc_OLc = [AfOc_OJc, AfOc_OUc];
    AfOc_OLc = sparse(AfOc_OLc);

    AOOLl = sparse([AfOl, AfOl_OLl;...
                zeros(size(AOLl,1), size(AfOl,1)), AOLl]);
    AOOLc = sparse([AfOc,  AfOc_OLc;...
                zeros(size(AOLc,1), size(AfOc,1)), AOLc]);
    AOOLl_OOLc = [param.pf.pcs.*speye(gridSf.a.n,gridSf.a.n),AfOl_OLc;...
            sparse(size(AOLl,1),gridSf.a.n),param.pcp.*speye(size(AOLl))];
    AOOL = [AOOLl, AOOLl_OOLc;...
               sparse(size(AOOLl_OOLc,1),size(AOOLl_OOLc,2)), AOOLc];
    solutionSP.AOOL = AOOL;     

    u_OOL = [uOl_f;uOLl;uOc_f;uOLc];


    %%===================================================================
    % Solve problem for the labor force Male nonparticipant
    %%===================================================================
    
    [AmOl_Ol, AmOl_Jl, AmOl_JOl, AmOl_UOl, AmOl_JOc, AmOl_UOc] = generateAOl_m(solutionSP, solutionSf_s, param.pm, gridSm);    
    [AmOc_Oc, AmOc_Jc, AmOc_JOc, AmOc_UOc] = generateAOc_m(solutionSP, solutionSf_s, param.pm, gridSm);    AmOl = AmOl_Ol;
    AmOc = AmOc_Oc;    
    AmOl = AmOl - param.pm.pcs.*speye(gridSm.a.n,gridSm.a.n);
    % AmOl = AmOl_Ol ...
    %      - param.pm.pcs .* speye(nOm) ...
    %      - spdiags(sum(AmOl_Jl,2), 0, nOm, nOm);      % Om -> Jm outflow
    % 
    % AmOc = AmOc_Oc ...
    %      - spdiags(sum(AmOc_Jc,2), 0, nOm, nOm);      % Om^c -> Jm^c outflow

    [AJOl_JOl, AJOl_UOl, AJOl_OOl, AJOl_JJl] = generateAJOl(solutionSP, param, gridP);
    [AUOl_UOl, AUOl_JOl, AUOl_OOl, AUOl_UJl] = generateAUOl(solutionSP, param, gridP);    
    [AJOc_JOc, AJOc_UOc, AJOc_OOc, AJOc_JJc] = generateAJOc(solutionSP, param, gridP);
    [AUOc_UOc, AUOc_JOc, AUOc_OOc, AUOc_UJc] = generateAUOc(solutionSP, param, gridP);    
        
    ALOl = [AJOl_JOl, AJOl_UOl;...
         AUOl_JOl, AUOl_UOl];
    ALOl = ALOl - param.pcp.*speye(param.NP/gridP.wf.n+gridP.af.n*gridP.am.n,param.NP/gridP.wf.n+gridP.af.n*gridP.am.n);
    ALOc = [AJOc_JOc, AJOc_UOc;...
         AUOc_JOc, AUOc_UOc];
    % ALOl = [AJOl_JOl, AJOl_UOl;
    %     AUOl_JOl, AUOl_UOl] ...
    %  - param.pcp .* speye(nJO+nUO);
    % ALOl = ALOl ...
    %  - blkdiag( spdiags(sum(AJOl_JJl,2),0,nJO,nJO), ... % JO_l -> JJ_l
    %             spdiags(sum(AUOl_UJl,2),0,nUO,nUO) );    % UO_l -> UJ_l
    % ALOc = [AJOc_JOc, AJOc_UOc;
    %         AUOc_JOc, AUOc_UOc];
    % ALOc = ALOc ...
    %      - blkdiag( spdiags(sum(AJOc_JJc,2),0,nJO,nJO), ... % JO_c -> JJ_c
    %                 spdiags(sum(AUOc_UJc,2),0,nUO,nUO) );    % UO_c -> UJ_c

    AmOl_LOl = [AmOl_JOl, AmOl_UOl];
    AmOl_LOl = sparse(AmOl_LOl);
    AmOl_LOc = [AmOl_JOc, AmOl_UOc];
    AmOl_LOc = sparse(AmOl_LOc);
    AmOc_LOc = [AmOc_JOc, AmOc_UOc];
    AmOc_LOc = sparse(AmOc_LOc);

    AOLOl = sparse([AmOl, AmOl_LOl;...
                zeros(size(ALOl,1), size(AmOl,1)), ALOl]);
    AOLOc = sparse([AmOc,  AmOc_LOc;...
                zeros(size(ALOc,1), size(AmOc,1)), ALOc]);
    AOLOl_OLOc = [param.pm.pcs.*speye(gridSm.a.n,gridSm.a.n),AmOl_LOc;...
            sparse(size(ALOl,1),gridSm.a.n),param.pcp.*speye(size(ALOl))];
    AOLO = [AOLOl, AOLOl_OLOc;...
               sparse(size(AOLOl_OLOc,1),size(AOLOl_OLOc,2)), AOLOc];
    solutionSP.AOLO = AOLO; 

    
    u_OLO = [uOl_m;uLOl;uOc_m;uLOc];


    %%===================================================================
    % Solve problem for the labor force Female and Male nonparticipant
    %%===================================================================
         
    [AOOl_OOl, AOOl_JOl, AOOl_OJl] = generateAOOl(solutionSP, param, gridP);
    [AOOc_OOc, AOOc_JOc, AOOc_OJc] = generateAOOc(solutionSP, param, gridP);
        
    AOO = [AOOl_OOl- param.pcp.*speye(gridP.af.n*gridP.am.n,gridP.af.n*gridP.am.n), param.pcp.*speye(gridP.af.n*gridP.am.n,gridP.af.n*gridP.am.n);...
         sparse(gridP.af.n*gridP.am.n,gridP.af.n*gridP.am.n), AOOc_OOc];
    % % ---- BOTH O (OO; no kids & kids) ----
    % D_OOl = spdiags(sum(AOOl_OJl,2) + sum(AOOl_JOl,2), 0, nOO, nOO);
    % D_OOc = spdiags(sum(AOOc_OJc,2) + sum(AOOc_JOc,2), 0, nOO, nOO);
    % AOO = [ (AOOl_OOl - D_OOl - param.pcp.*speye(nOO)),  param.pcp.*speye(nOO);
    %         sparse(nOO,nOO),                             (AOOc_OOc - D_OOc)     ];

    solutionSP.AOO = AOO;
       
    u_OO=[uOOl;uOOc];

%%===============================
% Build the off–diagonal transitions (index-safe)
%%===============================

% -------- sizes per half (no kids) --------
nLL_l  = (nJf+nUf+nJm+nUm+nJJ+nJU+nUJ+nUU);
nOOL_l = (nOf+nOJ+nOU);
nOLO_l = (nOm+nJO+nUO);
nOO_l  = (nOO);

% -------- global column starts (1-based) --------
C.LLL.l = 1;
C.LLL.c = C.LLL.l + nLL_l;

C.OOL.l = nALLL + 1;
C.OOL.c = C.OOL.l + nOOL_l;

C.OLO.l = nALLL + nAOOL + 1;
C.OLO.c = C.OLO.l + nOLO_l;

C.OO.l  = nALLL + nAOOL + nAOLO + 1;
C.OO.c  = C.OO.l + nOO_l;

% -------- LLL column sub-blocks --------
C.LLL.Jf = C.LLL.l;
C.LLL.Uf = C.LLL.Jf + nJf;
C.LLL.Jm = C.LLL.Uf + nUf;
C.LLL.Um = C.LLL.Jm + nJm;
C.LLL.JJ = C.LLL.Um + nUm;
C.LLL.JU = C.LLL.JJ + nJJ;
C.LLL.UJ = C.LLL.JU + nJU;
C.LLL.UU = C.LLL.UJ + nUJ;

C.LLLc.Jf = C.LLL.c;
C.LLLc.Uf = C.LLLc.Jf + nJf;
C.LLLc.Jm = C.LLLc.Uf + nUf;
C.LLLc.Um = C.LLLc.Jm + nJm;
C.LLLc.JJ = C.LLLc.Um + nUm;
C.LLLc.JU = C.LLLc.JJ + nJJ;
C.LLLc.UJ = C.LLLc.JU + nJU;
C.LLLc.UU = C.LLLc.UJ + nUJ;

% -------- OOL column sub-blocks --------
C.OOL.Of = C.OOL.l;
C.OOL.OJ = C.OOL.Of + nOf;
C.OOL.OU = C.OOL.OJ + nOJ;

C.OOLc.Of = C.OOL.c;
C.OOLc.OJ = C.OOLc.Of + nOf;
C.OOLc.OU = C.OOLc.OJ + nOJ;

% -------- OLO column sub-blocks --------
C.OLO.Om = C.OLO.l;
C.OLO.JO = C.OLO.Om + nOm;
C.OLO.UO = C.OLO.JO + nJO;

C.OLOc.Om = C.OLO.c;
C.OLOc.JO = C.OLOc.Om + nOm;
C.OLOc.UO = C.OLOc.JO + nJO;

% -------- OO column sub-blocks --------
C.OO.OOl = C.OO.l;
C.OO.OOc = C.OO.c;

% -------- row starts for each line (1-based) --------
R.LLL.Jf = 1;
R.LLL.Uf = R.LLL.Jf + nJf;
R.LLL.Jm = R.LLL.Uf + nUf;
R.LLL.Um = R.LLL.Jm + nJm;
R.LLL.JJ = R.LLL.Um + nUm;
R.LLL.JU = R.LLL.JJ + nJJ;
R.LLL.UJ = R.LLL.JU + nJU;
R.LLL.UU = R.LLL.UJ + nUJ;

% kids half of LLL starts exactly after the first half
R.LLLc.Jf = nLL_l + 1;
R.LLLc.Uf = R.LLLc.Jf + nJf;
R.LLLc.Jm = R.LLLc.Uf + nUf;
R.LLLc.Um = R.LLLc.Jm + nJm;
R.LLLc.JJ = R.LLLc.Um + nUm;
R.LLLc.JU = R.LLLc.JJ + nJJ;
R.LLLc.UJ = R.LLLc.JU + nJU;
R.LLLc.UU = R.LLLc.UJ + nUJ;

R.OOL.Of  = 1;
R.OOL.OJ  = R.OOL.Of + nOf;
R.OOL.OU  = R.OOL.OJ + nOJ;
R.OOLc.Of = nOOL_l + 1;
R.OOLc.OJ = R.OOLc.Of + nOf;
R.OOLc.OU = R.OOLc.OJ + nOJ;

R.OLO.Om  = 1;
R.OLO.JO  = R.OLO.Om + nOm;
R.OLO.UO  = R.OLO.JO + nJO;
R.OLOc.Om = nOLO_l + 1;
R.OLOc.JO = R.OLOc.Om + nOm;
R.OLOc.UO = R.OLOc.JO + nJO;

R.OO.OOl  = 1;
R.OO.OOc  = nOO_l + 1;

% Small range helper
irange = @(s, n) (s:(s+n-1));%(n>0) .* 

nTotCols = nALLL + nAOOL + nAOLO + nAOO;

%%========================
% LLL line
%%=======================
ALLLLine = sparse(nALLL, nTotCols);

% Place the big LLL block itself (into LLL columns)
ALLLLine(:, irange(C.LLL.l, nALLL)) = ALLL;

% no-kid single female → OOL / OLO
ALLLLine(irange(R.LLL.Jf, nJf), irange(C.OOL.Of, nOf)) = AfJl_Ol;
ALLLLine(irange(R.LLL.Jf, nJf), irange(C.OLO.JO, nJO)) = AfJl_JOl;
ALLLLine(irange(R.LLL.Jf, nJf), irange(C.OLOc.JO, nJO)) = AfJl_JOc;

ALLLLine(irange(R.LLL.Uf, nUf), irange(C.OOL.Of, nOf)) = AfUl_Ol;
ALLLLine(irange(R.LLL.Uf, nUf), irange(C.OLO.UO, nUO)) = AfUl_UOl;
ALLLLine(irange(R.LLL.Uf, nUf), irange(C.OLOc.UO, nUO)) = AfUl_UOc;

% no-kid single male → OOL / OLO
ALLLLine(irange(R.LLL.Jm, nJm), irange(C.OLO.Om, nOm)) = AmJl_Ol;   % Jm → Om(ℓ)
ALLLLine(irange(R.LLL.Jm, nJm), irange(C.OOL.OJ, nOJ)) = AmJl_OJl;  % Jm → OJ(ℓ)
ALLLLine(irange(R.LLL.Jm, nJm), irange(C.OOLc.OJ, nOJ)) = AmJl_OJc; % Jm → OJ(c)

ALLLLine(irange(R.LLL.Um, nUm), irange(C.OLO.Om, nOm)) = AmUl_Ol;   % Um → Om(ℓ)
ALLLLine(irange(R.LLL.Um, nUm), irange(C.OOL.OU, nOU)) = AmUl_OUl;  % Um → OU(ℓ)
ALLLLine(irange(R.LLL.Um, nUm), irange(C.OOLc.OU, nOU)) = AmUl_OUc; % Um → OU(c)

% no-kid couple → OOL / OLO
v0 = R.LLL.JJ;
ALLLLine(irange(v0, nJJ), irange(C.OOL.OJ, nOJ)) = AJJl_OJl;
ALLLLine(irange(v0, nJJ), irange(C.OLO.JO, nJO)) = AJJl_JOl;

v0 = R.LLL.JU;
ALLLLine(irange(v0, nJU), irange(C.OOL.OU, nOU)) = AJUl_OUl;
ALLLLine(irange(v0, nJU), irange(C.OLO.JO, nJO)) = AJUl_JOl;

v0 = R.LLL.UJ;
ALLLLine(irange(v0, nUJ), irange(C.OOL.OJ, nOJ)) = AUJl_OJl;
ALLLLine(irange(v0, nUJ), irange(C.OLO.UO, nUO)) = AUJl_UOl;

v0 = R.LLL.UU;
ALLLLine(irange(v0, nUU), irange(C.OOL.OU, nOU)) = AUUl_OUl;
ALLLLine(irange(v0, nUU), irange(C.OLO.UO, nUO)) = AUUl_UOl;

% kid single female/male → OOLc / OLOc
v0 = R.LLLc.Jf;
ALLLLine(irange(v0, nJf), irange(C.OOLc.Of, nOf)) = AfJc_Oc;
ALLLLine(irange(v0, nJf), irange(C.OLOc.JO, nJO)) = AfJc_JOc;

v0 = R.LLLc.Uf;
ALLLLine(irange(v0, nUf), irange(C.OOLc.Of, nOf)) = AfUc_Oc;
ALLLLine(irange(v0, nUf), irange(C.OLOc.UO, nUO)) = AfUc_UOc;

v0 = R.LLLc.Jm;
ALLLLine(irange(v0, nJm), irange(C.OLOc.Om, nOm)) = AmJc_Oc;
ALLLLine(irange(v0, nJm), irange(C.OOLc.OJ, nOJ)) = AmJc_OJc;

v0 = R.LLLc.Um;
ALLLLine(irange(v0, nUm), irange(C.OLOc.Om, nOm)) = AmUc_Oc;
ALLLLine(irange(v0, nUm), irange(C.OOLc.OU, nOU)) = AmUc_OUc;

% kid couple → OOLc / OLOc
v0 = R.LLLc.JJ;
ALLLLine(irange(v0, nJJ), irange(C.OOLc.OJ, nOJ)) = AJJc_OJc;
ALLLLine(irange(v0, nJJ), irange(C.OLOc.JO, nJO)) = AJJc_JOc;

v0 = R.LLLc.JU;
ALLLLine(irange(v0, nJU), irange(C.OOLc.OU, nOU)) = AJUc_OUc;
ALLLLine(irange(v0, nJU), irange(C.OLOc.JO, nJO)) = AJUc_JOc;

v0 = R.LLLc.UJ;
ALLLLine(irange(v0, nUJ), irange(C.OOLc.OJ, nOJ)) = AUJc_OJc;
ALLLLine(irange(v0, nUJ), irange(C.OLOc.UO, nUO)) = AUJc_UOc;

v0 = R.LLLc.UU;
ALLLLine(irange(v0, nUU), irange(C.OOLc.OU, nOU)) = AUUc_OUc;
ALLLLine(irange(v0, nUU), irange(C.OLOc.UO, nUO)) = AUUc_UOc;

%%========================
% OOL line (female O)
%%=======================
AOOLLine = sparse(nAOOL, nTotCols);

% OOL on its own columns
AOOLLine(:, irange(C.OOL.l, nAOOL)) = AOOL;

% no-kid: Of → Jf ; OJ → JJ ; OU → JU
AOOLLine(irange(R.OOL.Of, nOf), irange(C.LLL.Jf, nJf)) = AfOl_Jl;
AOOLLine(irange(R.OOL.OJ, nOJ), irange(C.LLL.JJ, nJJ)) = AOJl_JJl;
AOOLLine(irange(R.OOL.OU, nOU), irange(C.LLL.JU, nJU)) = AOUl_JUl;

% no-kid to OO
AOOLLine(irange(R.OOL.OJ, nOJ), irange(C.OO.OOl, nOO)) = AOJl_OOl;
AOOLLine(irange(R.OOL.OU, nOU), irange(C.OO.OOl, nOO)) = AOUl_OOl;

% kids: Ofᶜ → Jfᶜ ; OJᶜ → JJᶜ ; OUᶜ → JUᶜ
AOOLLine(irange(R.OOLc.Of, nOf), irange(C.LLLc.Jf, nJf)) = AfOc_Jc;
AOOLLine(irange(R.OOLc.OJ, nOJ), irange(C.LLLc.JJ, nJJ)) = AOJc_JJc;
AOOLLine(irange(R.OOLc.OU, nOU), irange(C.LLLc.JU, nJU)) = AOUc_JUc;

% kids to OOᶜ
AOOLLine(irange(R.OOLc.OJ, nOJ), irange(C.OO.OOc, nOO)) = AOJc_OOc;
AOOLLine(irange(R.OOLc.OU, nOU), irange(C.OO.OOc, nOO)) = AOUc_OOc;

%%========================
% OLO line (male O)
%%=======================
AOLOLine = sparse(nAOLO, nTotCols);

% OLO on its own columns
AOLOLine(:, irange(C.OLO.l, nAOLO)) = AOLO;

% no-kid: Om → Jm ; JO → JJ ; UO → UJ
AOLOLine(irange(R.OLO.Om, nOm), irange(C.LLL.Jm, nJm)) = AmOl_Jl;
AOLOLine(irange(R.OLO.JO, nJO), irange(C.LLL.JJ, nJJ)) = AJOl_JJl;
AOLOLine(irange(R.OLO.UO, nUO), irange(C.LLL.UJ, nUJ)) = AUOl_UJl;

% no-kid to OO
AOLOLine(irange(R.OLO.JO, nJO), irange(C.OO.OOl, nOO)) = AJOl_OOl;
AOLOLine(irange(R.OLO.UO, nUO), irange(C.OO.OOl, nOO)) = AUOl_OOl;

% kids: Omᶜ → Jmᶜ ; JOᶜ → JJᶜ ; UOᶜ → UJᶜ
AOLOLine(irange(R.OLOc.Om, nOm), irange(C.LLLc.Jm, nJm)) = AmOc_Jc;
AOLOLine(irange(R.OLOc.JO, nJO), irange(C.LLLc.JJ, nJJ)) = AJOc_JJc;
AOLOLine(irange(R.OLOc.UO, nUO), irange(C.LLLc.UJ, nUJ)) = AUOc_UJc;

% kids to OOᶜ
AOLOLine(irange(R.OLOc.JO, nJO), irange(C.OO.OOc, nOO)) = AJOc_OOc;
AOLOLine(irange(R.OLOc.UO, nUO), irange(C.OO.OOc, nOO)) = AUOc_OOc;

%%========================
% OO line (both O)
%%=======================
AOOLine = sparse(nAOO, nTotCols);

% OO on its own columns
AOOLine(:, irange(C.OO.l, nAOO)) = AOO;

% no-kid OO → OJ(ℓ)/JO(ℓ)
AOOLine(irange(R.OO.OOl, nOO), irange(C.OOL.OJ, nOJ)) = AOOl_OJl;
AOOLine(irange(R.OO.OOl, nOO), irange(C.OLO.JO, nJO)) = AOOl_JOl;

% kids OOᶜ → OJᶜ / UJOᶜ
AOOLine(irange(R.OO.OOc, nOO), irange(C.OOLc.OJ, nOJ)) = AOOc_OJc;
AOOLine(irange(R.OO.OOc, nOO), irange(C.OLOc.JO, nJO)) = AOOc_JOc;

%%========================
% Stack lines
%%=======================
A = [ALLLLine;
     AOOLLine;
     AOLOLine;
     AOOLine];

solutionSP.Aall = A;

solutionSP.Aall = A;
[nr, nc] = size(A);
assert(nr == nc && nr == nTotCols, 'A must be square with nTotCols columns.');

% --------- HARDEN A AS A CONTINUOUS-TIME GENERATOR ----------
tol = 1e-12;

% (a) Kill tiny negative off-diagonals (keep the diagonal separate)
dA   = spdiags(A,0);
Aoff = A - spdiags(dA,0,nr,nc);
negMask = Aoff < -tol;
if nnz(negMask)
    % Clamp only the genuinely tiny negatives
    Aoff(negMask) = 0;
end
A = Aoff + spdiags(dA,0,nr,nc);

% (b) Force exact row-sum = 0  -> diag(A) = -sum(offdiag)
rowSum = sum(A,2);                 % includes diagonal
A = A - spdiags(rowSum, 0, nr, nc);

% (c) Sanity: diag(A) ≤ 0 and off-diagonals ≥ 0
if any(full(diag(A)) > tol)
    error('A has positive diagonal entries after correction.');
end
Aoff = A - spdiags(spdiags(A,0),0,nr,nc);
if any(Aoff(:) < -tol)
    error('A has negative off-diagonals after correction.');
end

% (d) Discounting. Step (b) zeroes every row sum, which also removes the
% retirement exit (-zeta) that the C-blocks put on the diagonal. Retirement
% is an exit to a zero-value absorbing state, so it is restored here as
% extra discounting: (rho + zeta) V = u + A V, as in the paper's HJB.
rho_use = param.rho + param.zeta;

B = rho_use .* speye(nr) - A;
u = [u_LLL; u_OOL; u_OLO; u_OO];

% (optional) quick condition diagnostics
% fprintf('condest(B) ~ %.2e\n', condest(B));

Vnew = B \ u;

% Optional quick sanity check:
% rs = full(sum(A,2)); fprintf('Row-sum min=%g, max=%g\n', min(rs), max(rs));

    % ---------- unpack values
    l=1; r=size(ALLL,1);     LLLnew = Vnew(l:r);
    l=r+1; r=r+size(AOOL,1); OOLnew = Vnew(l:r);
    l=r+1; r=r+size(AOLO,1); OLOnew = Vnew(l:r);
    l=r+1; r=r+size(AOO,1);  OOnew  = Vnew(l:r);

    solutionSP.L_LLvec = LLLnew;

    % LLL unpack
    high = 0;
    low=1; high=high+nJf; JflValueNew = reshape(LLLnew(low:high), size(solutionSP.JflValue));
    low=high+1; high=high+nUf; UflValueNew = LLLnew(low:high);
    low=high+1; high=high+nJm; JmlValueNew = reshape(LLLnew(low:high), size(solutionSP.JmlValue));
    low=high+1; high=high+nUm; UmValueNew  = LLLnew(low:high);

    low=high+1; high=high+nJJ; JJlValueNew = reshape(LLLnew(low: high), size(solutionSP.JJlValue));
    low=high+1; high=high+nJU; JUlValueNew = reshape(LLLnew(low:high), size(solutionSP.JUlValue));
    low=high+1; high=high+nUJ; UJlValueNew = reshape(LLLnew(low:high), size(solutionSP.UJlValue));
    low=high+1; high=high+nUU; UUlValueNew = reshape(LLLnew(low:high), size(solutionSP.UUlValue));

    m = high;
    low=m+1; high=high+nJf; JfcValueNew = reshape(LLLnew(low:high), size(solutionSP.JfcValue));
    low=high+1; high=high+nUf; UfcValueNew = LLLnew(low:high);
    low=high+1; high=high+nJm; JmcValueNew = reshape(LLLnew(low:high), size(solutionSP.JmcValue));
    low=high+1; high=high+nUm; UmcValueNew = LLLnew(low:high);

    low=high+1; high=high+nJJ; JJcValueNew = reshape(LLLnew(low: high), size(solutionSP.JJcValue));
    low=high+1; high=high+nJU; JUcValueNew = reshape(LLLnew(low:high), size(solutionSP.JUcValue));
    low=high+1; high=high+nUJ; UJcValueNew = reshape(LLLnew(low:high), size(solutionSP.UJcValue));
    low=high+1; high=high+nUU; UUcValueNew = reshape(LLLnew(low:high), size(solutionSP.UUcValue));

    diffJfl = max(abs(JflValueNew-solutionSP.JflValue),[],'all');
    diffUfl = max(abs(UflValueNew-solutionSP.UflValue),[],'all');
    diffJml = max(abs(JmlValueNew-solutionSP.JmlValue),[],'all');
    diffUml = max(abs(UmValueNew -solutionSP.UmlValue),[],'all');
    diffJJl = max(abs(JJlValueNew-solutionSP.JJlValue),[],'all');
    diffJUl = max(abs(JUlValueNew-solutionSP.JUlValue),[],'all');
    diffUJl = max(abs(UJlValueNew-solutionSP.UJlValue),[],'all');
    diffUUl = max(abs(UUlValueNew-solutionSP.UUlValue),[],'all');

    diffJfc = max(abs(JfcValueNew-solutionSP.JfcValue),[],'all');
    diffUfc = max(abs(UfcValueNew-solutionSP.UfcValue),[],'all');
    diffJmc = max(abs(JmcValueNew-solutionSP.JmcValue),[],'all');
    diffUmc = max(abs(UmcValueNew-solutionSP.UmcValue),[],'all');
    diffJJc = max(abs(JJcValueNew-solutionSP.JJcValue),[],'all');
    diffJUc = max(abs(JUcValueNew-solutionSP.JUcValue),[],'all');
    diffUJc = max(abs(UJcValueNew-solutionSP.UJcValue),[],'all');
    diffUUc = max(abs(UUcValueNew-solutionSP.UUcValue),[],'all');

    diffL  = max([diffJfl, diffUfl, diffJfc, diffUfc, diffJml, diffUml, diffJmc, diffUmc]);
    diffPLL= max([diffJJl, diffJUl, diffUJl, diffUUl, diffJJc, diffJUc, diffUJc, diffUUc]);

    solutionSP.JflValue = (1-param.updMult_s ).*solutionSP.JflValue + param.updMult_s.*JflValueNew;
    solutionSP.UflValue = (1-param.updMult_s ).*solutionSP.UflValue + param.updMult_s.*UflValueNew;
    solutionSP.JmlValue = (1-param.updMult_s ).*solutionSP.JmlValue + param.updMult_s.*JmlValueNew;
    solutionSP.UmlValue = (1-param.updMult_s ).*solutionSP.UmlValue + param.updMult_s.*UmValueNew;

    solutionSP.JJlValue = (1-param.updMult_d ).*solutionSP.JJlValue + param.updMult_d.*JJlValueNew;
    solutionSP.JUlValue = (1-param.updMult_d ).*solutionSP.JUlValue + param.updMult_d.*JUlValueNew;
    solutionSP.UJlValue = (1-param.updMult_d ).*solutionSP.UJlValue + param.updMult_d.*UJlValueNew;
    solutionSP.UUlValue = (1-param.updMult_d ).*solutionSP.UUlValue + param.updMult_d.*UUlValueNew;

    solutionSP.JfcValue = (1-param.updMult_s ).*solutionSP.JfcValue + param.updMult_s.*JfcValueNew;
    solutionSP.UfcValue = (1-param.updMult_s ).*solutionSP.UfcValue + param.updMult_s.*UfcValueNew;
    solutionSP.JmcValue = (1-param.updMult_s ).*solutionSP.JmcValue + param.updMult_s.*JmcValueNew;
    solutionSP.UmcValue = (1-param.updMult_s ).*solutionSP.UmcValue + param.updMult_s.*UmcValueNew;

    solutionSP.JJcValue = (1-param.updMult_d ).*solutionSP.JJcValue + param.updMult_d.*JJcValueNew;
    solutionSP.JUcValue = (1-param.updMult_d ).*solutionSP.JUcValue + param.updMult_d.*JUcValueNew;
    solutionSP.UJcValue = (1-param.updMult_d ).*solutionSP.UJcValue + param.updMult_d.*UJcValueNew;
    solutionSP.UUcValue = (1-param.updMult_d ).*solutionSP.UUcValue + param.updMult_d.*UUcValueNew;

    % OOL unpack
    solutionSP.O_OLvec = OOLnew;
    high=0;
    low=1; high=high+nOf; OflValueNew = OOLnew(low:high);
    low=high+1; high=high+nOJ; OJlValueNew = reshape(OOLnew(low:high), size(solutionSP.OJlValue));
    low=high+1; high=high+nOU; OUlValueNew = reshape(OOLnew(low:high), size(solutionSP.OUlValue));

    m=high;
    low=m+1; high=high+nOf; OfcValueNew = OOLnew(low:high);
    low=high+1; high=high+nOJ; OJcValueNew = reshape(OOLnew(low:high), size(solutionSP.OJcValue));
    low=high+1; high=high+nOU; OUcValueNew = reshape(OOLnew(low:high), size(solutionSP.OUcValue));

    diffOf  = max(abs(OflValueNew-solutionSP.OflValue),[],'all');
    diffOJ  = max(abs(OJlValueNew-solutionSP.OJlValue),[],'all');
    diffOU  = max(abs(OUlValueNew-solutionSP.OUlValue),[],'all');
    diffOfc = max(abs(OfcValueNew-solutionSP.OfcValue),[],'all');
    diffOJc = max(abs(OJcValueNew-solutionSP.OJcValue),[],'all');
    diffOUc = max(abs(OUcValueNew-solutionSP.OUcValue),[],'all');

    diffSOf = max([diffOf,diffOfc]);
    diffPOL = max([diffOJ, diffOU, diffOJc, diffOUc]);

    solutionSP.OflValue = (1-param.updMult_s ).*solutionSP.OflValue + param.updMult_s.*OflValueNew;
    solutionSP.OJlValue = (1-param.updMult_d ).*solutionSP.OJlValue + param.updMult_d.*OJlValueNew;
    solutionSP.OUlValue = (1-param.updMult_d ).*solutionSP.OUlValue + param.updMult_d.*OUlValueNew;
    solutionSP.OfcValue = (1-param.updMult_s ).*solutionSP.OfcValue + param.updMult_s.*OfcValueNew;
    solutionSP.OJcValue = (1-param.updMult_d ).*solutionSP.OJcValue + param.updMult_d.*OJcValueNew;
    solutionSP.OUcValue = (1-param.updMult_d ).*solutionSP.OUcValue + param.updMult_d.*OUcValueNew;

    % OLO unpack
    solutionSP.O_LOvec = OLOnew;
    high=0;
    low=1; high=high+nOm; OmlValueNew = OLOnew(low:high);
    low=high+1; high=high+nJO; JOlValueNew = reshape(OLOnew(low:high), size(solutionSP.JOlValue));
    low=high+1; high=high+nUO; UOlValueNew = reshape(OLOnew(low:high), size(solutionSP.UOlValue));

    m=high;
    low=m+1; high=high+nOm; OmcValueNew = OLOnew(low:high);
    low=high+1; high=high+nJO; JOcValueNew = reshape(OLOnew(low:high), size(solutionSP.JOcValue));
    low=high+1; high=high+nUO; UOcValueNew = reshape(OLOnew(low:high), size(solutionSP.UOcValue));

    diffOml = max(abs(OmlValueNew-solutionSP.OmlValue),[],'all');
    diffJOl = max(abs(JOlValueNew-solutionSP.JOlValue),[],'all');
    diffUOl = max(abs(UOlValueNew-solutionSP.UOlValue),[],'all');
    diffOmc = max(abs(OmcValueNew-solutionSP.OmcValue),[],'all');
    diffJOc = max(abs(JOcValueNew-solutionSP.JOcValue),[],'all');
    diffUOc = max(abs(UOcValueNew-solutionSP.UOcValue),[],'all');

    diffSOm = max([diffOml,diffOmc]);
    diffPLO = max([diffJOl, diffUOl, diffJOc, diffUOc]);

    solutionSP.OmlValue = (1-param.updMult_s ).*solutionSP.OmlValue + param.updMult_s.*OmlValueNew;
    solutionSP.JOlValue = (1-param.updMult_d ).*solutionSP.JOlValue + param.updMult_d.*JOlValueNew;
    solutionSP.UOlValue = (1-param.updMult_d ).*solutionSP.UOlValue + param.updMult_d.*UOlValueNew;
    solutionSP.OmcValue = (1-param.updMult_s ).*solutionSP.OmcValue + param.updMult_s.*OmcValueNew;
    solutionSP.JOcValue = (1-param.updMult_d ).*solutionSP.JOcValue + param.updMult_d.*JOcValueNew;
    solutionSP.UOcValue = (1-param.updMult_d ).*solutionSP.UOcValue + param.updMult_d.*UOcValueNew;

    % OO unpack
    solutionSP.OOvec = OOnew;
    OOlValueNew = reshape(OOnew(1:nOO), size(solutionSP.OOlValue));
    OOcValueNew = reshape(OOnew(nOO+1:end), size(solutionSP.OOcValue));

    diffOOl = max(abs(OOlValueNew-solutionSP.OOlValue),[],'all');
    diffOOc = max(abs(OOcValueNew-solutionSP.OOcValue),[],'all');
    diffPOO  = max([diffOOl, diffOOc]);

    solutionSP.OOlValue = (1-param.updMult_d).*solutionSP.OOlValue+param.updMult_d.*OOlValueNew;
    solutionSP.OOcValue = (1-param.updMult_d).*solutionSP.OOcValue+param.updMult_d.*OOcValueNew;  
    
end


% functions for single
function [AJ_J, AJ_U, AJ_O, AJ_JJ, AJ_JU, AJ_JO, AJ_JJc, AJ_JUc, AJ_JOc] = generateAJl_f(solutionSP,solutionS_s, pg, gridS)
    AJ_J = (zeros(pg.NS,pg.NS));
    AJ_O = (zeros(size(pg.CJ_U)));
    AJ_JJ = (zeros(pg.NS,param.NP));
    AJ_JU = (zeros(pg.NS,param.NP/gridP.wm.n));
    AJ_JO = (zeros(pg.NS,param.NP/gridP.wm.n));
    AJ_JJc = (zeros(pg.NS,param.NP));
    AJ_JUc = (zeros(pg.NS,param.NP/gridP.wm.n));
    AJ_JOc = (zeros(pg.NS,param.NP/gridP.wm.n));
    singles=1;%sum(solutionSP.h,'all')+sum(solutionSP.e,'all');

    for ai = 1:gridS.a.n
        %find weights for linear interpolation
        ai_m=fVai_m(ai,:);
        ai_m_w=fVai_m_w(ai,:);

        for wi = 1:gridS.w.n
            %find weights for linear interpolation
            wi_m=fVwi_m(wi,:);
            wi_m_w=fVwi_m_w(wi,:);

            i = (wi-1)*gridS.a.n + ai;
            if solutionSP.JflValue(ai, wi) < solutionSP.OflValue(ai)
                AJ_J(i,i)   = AJ_J(i,i)   - pg.lambda_o;
                AJ_O(i,ai)  = AJ_O(i,ai)  + pg.lambda_o;
            elseif tremble~=0
                t = tremble/exp( solutionSP.JflValue(ai, wi) - solutionSP.OflValue(ai) );
                AJ_J(i,i)   = AJ_J(i,i)   - t;
                AJ_O(i,ai)  = AJ_O(i,ai)  + t;
            end

            % for wi_new = 1:gridS.w.n
            %     i_new = (wi_new-1)*gridS.a.n + ai;
            %     if solutionSP.JflValue(ai, wi)<solutionSP.JflValue(ai, wi_new)
            %         AJ_J(i,i) = AJ_J(i,i)-pg.lambda*param.f_long(wi_new);
            %         AJ_J(i,i_new) = AJ_J(i,i_new)+pg.lambda*param.f_long(wi_new);generateAJl_f
            %     end
            % end
            for wi_new = 1:gridP.wm.n
                %====
                %check marrige
                for ai_new = 1:gridP.am.n                      
                    JJValue=0;                    
                    JJcValue=0;                    
                    for k=1:4
                        sel='00';
                        binary = dec2bin(k-1);
                        sel(end-numel(binary)+1:end)=binary;
                        sel=reshape(str2num(char(num2cell(sel))),1,[]);
                        sel=sel+1;
                        JJValue=JJValue+...
                            solutionSP.JJlValue(ai_m(sel(1)), wi_m(sel(2)), ai_new, wi_new)*...
                            (ai_m_w(sel(1))*wi_m_w(sel(2)));                        
                        JJcValue=JJcValue+...
                            solutionSP.JJcValue(ai_m(sel(1)), wi_m(sel(2)), ai_new, wi_new)*...
                            (ai_m_w(sel(1))*wi_m_w(sel(2)));                        
                    end                              

                    if always_marry||((solutionSP.JflValue(ai, wi)<JJValue)&&...
                       (solutionSP.JmlValue((ai_new-1)*param.Sdimmult+1, (wi_new-1)*param.Sdimmult+1)<JJValue))
                        for k=1:4
                            sel='00';
                            binary = dec2bin(k-1);
                            sel(end-numel(binary)+1:end)=binary;
                            sel=reshape(str2num(char(num2cell(sel))),1,[]);
                            sel=sel+1;
                            i_self= (wi_new-1)*gridP.am.n*gridP.af.n*gridP.wf.n+...
                                (ai_new-1)*gridP.af.n*gridP.wf.n+(wi_m(sel(2))-1)*gridP.af.n+ai_m(sel(1));
                            AJ_J(i,i) = AJ_J(i,i)-param.m*solutionS_s.el(ai_new,wi_new)^1*(ai_m_w(sel(1))*wi_m_w(sel(2)))/singles;
                            AJ_JJ(i,i_self)=AJ_JJ(i,i_self)+param.m*solutionS_s.el(ai_new,wi_new)^1*(ai_m_w(sel(1))*wi_m_w(sel(2)))/singles;
                        end                        
                    end
                    if always_marry||((solutionSP.JflValue(ai, wi)<JJcValue)&&...
                       (solutionSP.JmcValue((ai_new-1)*param.Sdimmult+1, (wi_new-1)*param.Sdimmult+1)<JJcValue))
                        for k=1:4
                            sel='00';
                            binary = dec2bin(k-1);
                            sel(end-numel(binary)+1:end)=binary;
                            sel=reshape(str2num(char(num2cell(sel))),1,[]);
                            sel=sel+1;
                            i_self= (wi_new-1)*gridP.am.n*gridP.af.n*gridP.wf.n+...
                                (ai_new-1)*gridP.af.n*gridP.wf.n+(wi_m(sel(2))-1)*gridP.af.n+ai_m(sel(1));
                            AJ_J(i,i) = AJ_J(i,i)-param.m*solutionS_s.ec(ai_new,wi_new)^1*(ai_m_w(sel(1))*wi_m_w(sel(2)))/singles;
                            AJ_JJc(i,i_self)=AJ_JJc(i,i_self)+param.m*solutionS_s.ec(ai_new,wi_new)^1*(ai_m_w(sel(1))*wi_m_w(sel(2)))/singles;
                        end   
                    end
                    if wi_new==1

                        JUValue=0;
                        JUcValue=0;
                        for k=1:4
                            sel='00';
                            binary = dec2bin(k-1);
                            sel(end-numel(binary)+1:end)=binary;
                            sel=reshape(str2num(char(num2cell(sel))),1,[]);
                            sel=sel+1;
                            JUValue=JUValue+...
                                solutionSP.JUlValue(ai_m(sel(1)), wi_m(sel(2)), ai_new)*...
                                (ai_m_w(sel(1))*wi_m_w(sel(2)));                            
                            JUcValue=JUcValue+...
                                solutionSP.JUcValue(ai_m(sel(1)), wi_m(sel(2)), ai_new)*...
                                (ai_m_w(sel(1))*wi_m_w(sel(2)));                            
                        end


                        if always_marry||((solutionSP.JflValue(ai, wi)<JUValue)&&...
                           (solutionSP.UmlValue((ai_new-1)*param.Sdimmult+1)<JUValue))
                            for k=1:4
                                sel='00';
                                binary = dec2bin(k-1);
                                sel(end-numel(binary)+1:end)=binary;
                                sel=reshape(str2num(char(num2cell(sel))),1,[]);
                                sel=sel+1;                  
                                i_self= (ai_new-1)*gridP.af.n*gridP.wf.n+(wi_m(sel(2))-1)*gridP.af.n+ai_m(sel(1));
                                AJ_J(i,i) = AJ_J(i,i)-param.m*solutionS_s.hl(ai_new)^1*(ai_m_w(sel(1))*wi_m_w(sel(2)))/singles;
                                AJ_JU(i,i_self)=AJ_JU(i,i_self)+param.m*solutionS_s.hl(ai_new)^1*(ai_m_w(sel(1))*wi_m_w(sel(2)))/singles;
                            end  

                        end
                        if always_marry||((solutionSP.JflValue(ai, wi)<JUcValue)&&...
                           (solutionSP.UmcValue((ai_new-1)*param.Sdimmult+1)<JUcValue))
                            for k=1:4
                                sel='00';
                                binary = dec2bin(k-1);
                                sel(end-numel(binary)+1:end)=binary;
                                sel=reshape(str2num(char(num2cell(sel))),1,[]);
                                sel=sel+1;                  
                                i_self= (ai_new-1)*gridP.af.n*gridP.wf.n+(wi_m(sel(2))-1)*gridP.af.n+ai_m(sel(1));
                                AJ_J(i,i) = AJ_J(i,i)-param.m*solutionS_s.hc(ai_new)^1*(ai_m_w(sel(1))*wi_m_w(sel(2)))/singles;
                                AJ_JUc(i,i_self)=AJ_JUc(i,i_self)+param.m*solutionS_s.hc(ai_new)^1*(ai_m_w(sel(1))*wi_m_w(sel(2)))/singles;
                            end      
                        end

                        JOValue=0;
                        JOcValue=0;
                        for k=1:4
                            sel='00';
                            binary = dec2bin(k-1);
                            sel(end-numel(binary)+1:end)=binary;
                            sel=reshape(str2num(char(num2cell(sel))),1,[]);
                            sel=sel+1;
                            JOValue=JOValue+...
                                solutionSP.JOlValue(ai_m(sel(1)), wi_m(sel(2)), ai_new)*...
                                (ai_m_w(sel(1))*wi_m_w(sel(2)));                            
                            JOcValue=JOcValue+...
                                solutionSP.JOcValue(ai_m(sel(1)), wi_m(sel(2)), ai_new)*...
                                (ai_m_w(sel(1))*wi_m_w(sel(2)));                            
                        end
                        if always_marry||((solutionSP.JflValue(ai, wi)<JOValue)&&...
                           (solutionSP.OmlValue((ai_new-1)*param.Sdimmult+1)<JOValue))
                            for k=1:4
                                sel='00';
                                binary = dec2bin(k-1);
                                sel(end-numel(binary)+1:end)=binary;
                                sel=reshape(str2num(char(num2cell(sel))),1,[]);
                                sel=sel+1;                  
                                i_self= (ai_new-1)*gridP.af.n*gridP.wf.n+(wi_m(sel(2))-1)*gridP.af.n+ai_m(sel(1));
                                AJ_J(i,i) = AJ_J(i,i)-param.m*solutionS_s.ol(ai_new)^1*(ai_m_w(sel(1))*wi_m_w(sel(2)))/singles;
                                AJ_JO(i,i_self)=AJ_JO(i,i_self)+param.m*solutionS_s.ol(ai_new)^1*(ai_m_w(sel(1))*wi_m_w(sel(2)))/singles;
                            end  

                        end
                        if always_marry||((solutionSP.JflValue(ai, wi)<JOcValue)&&...
                           (solutionSP.OmcValue((ai_new-1)*param.Sdimmult+1)<JOcValue))
                            for k=1:4
                                sel='00';
                                binary = dec2bin(k-1);
                                sel(end-numel(binary)+1:end)=binary;
                                sel=reshape(str2num(char(num2cell(sel))),1,[]);
                                sel=sel+1;                  
                                i_self= (ai_new-1)*gridP.af.n*gridP.wf.n+(wi_m(sel(2))-1)*gridP.af.n+ai_m(sel(1));
                                AJ_J(i,i) = AJ_J(i,i)-param.m*solutionS_s.oc(ai_new)^1*(ai_m_w(sel(1))*wi_m_w(sel(2)))/singles;
                                AJ_JOc(i,i_self)=AJ_JOc(i,i_self)+param.m*solutionS_s.oc(ai_new)^1*(ai_m_w(sel(1))*wi_m_w(sel(2)))/singles;
                            end      
                        end

                    end                
                end
            end
        end
    end   
    AJ_J = sparse(AJ_J)+pg.CJ;
    AJ_U = pg.CJ_U;  
    AJ_O = sparse(AJ_O);  
    AJ_JJ = sparse(AJ_JJ);
    AJ_JU = sparse(AJ_JU);
    AJ_JO = sparse(AJ_JO);
    AJ_JJc = sparse(AJ_JJc);
    AJ_JUc = sparse(AJ_JUc);
    AJ_JOc = sparse(AJ_JOc);
end

function [AJ_J, AJ_U, AJ_O, AJ_JJ, AJ_UJ, AJ_OJ, AJ_JJc, AJ_UJc, AJ_OJc] = generateAJl_m(solutionSP,solutionS_s, pg, gridS)
    AJ_J = (zeros(pg.NS,pg.NS));
    AJ_O = zeros(size(pg.CJ_U));
    AJ_JJ = (zeros(pg.NS,param.NP));
    AJ_UJ = (zeros(pg.NS,param.NP/gridP.wf.n));
    AJ_OJ = (zeros(pg.NS,param.NP/gridP.wf.n));
    AJ_JJc = (zeros(pg.NS,param.NP));
    AJ_UJc = (zeros(pg.NS,param.NP/gridP.wf.n));
    AJ_OJc = (zeros(pg.NS,param.NP/gridP.wf.n));
    singles=1;%sum(solutionSP.h,'all')+sum(solutionSP.e,'all');

    for ai = 1:gridS.a.n
        %find weights for linear interpolation
        ai_m=mVai_m(ai,:);
        ai_m_w=mVai_m_w(ai,:);

        for wi = 1:gridS.w.n
            %find weights for linear interpolation
            wi_m=mVwi_m(wi,:);
            wi_m_w=mVwi_m_w(wi,:);

            i = (wi-1)*gridS.a.n + ai;
            if solutionSP.JmlValue(ai, wi) < solutionSP.OmlValue(ai)
                AJ_J(i,i)   = AJ_J(i,i)   - pg.lambda_o;
                AJ_O(i,ai)  = AJ_O(i,ai)  + pg.lambda_o;
            elseif tremble~=0
                t = tremble/exp( solutionSP.JmlValue(ai, wi) - solutionSP.OmlValue(ai) );
                AJ_J(i,i)   = AJ_J(i,i)   - t;
                AJ_O(i,ai)  = AJ_O(i,ai)  + t;
            end
            % for wi_new = 1:gridS.w.n
            %     i_new = (wi_new-1)*gridS.a.n + ai;
            %     if solutionSP.JmlValue(ai, wi)<solutionSP.JmlValue(ai, wi_new)
            %         AJ_J(i,i) = AJ_J(i,i)-pg.lambda*param.f_long(wi_new);
            %         AJ_J(i,i_new) = AJ_J(i,i_new)+pg.lambda*param.f_long(wi_new);
            %     end
            % end
            for wi_new = 1:gridP.wf.n
                %====
                %check marrige
                for ai_new = 1:gridP.af.n                      
                    JJValue=0;                    
                    JJcValue=0;                    
                    for k=1:4
                        sel='00';
                        binary = dec2bin(k-1);
                        sel(end-numel(binary)+1:end)=binary;
                        sel=reshape(str2num(char(num2cell(sel))),1,[]);
                        sel=sel+1;
                        JJValue=JJValue+...
                            solutionSP.JJlValue(ai_new, wi_new, ai_m(sel(1)), wi_m(sel(2)))*...
                            (ai_m_w(sel(1))*wi_m_w(sel(2)));                        
                        JJcValue=JJcValue+...
                            solutionSP.JJcValue(ai_new, wi_new, ai_m(sel(1)), wi_m(sel(2)))*...
                            (ai_m_w(sel(1))*wi_m_w(sel(2)));                        
                    end


                    if always_marry||((solutionSP.JmlValue(ai, wi)<JJValue)&&...
                       (solutionSP.JflValue((ai_new-1)*param.Sdimmult+1, (wi_new-1)*param.Sdimmult+1)<JJValue))
                        for k=1:4
                            sel='00';
                            binary = dec2bin(k-1);
                            sel(end-numel(binary)+1:end)=binary;
                            sel=reshape(str2num(char(num2cell(sel))),1,[]);
                            sel=sel+1;
                            i_self= (wi_m(sel(2))-1)*gridP.am.n*gridP.af.n*gridP.wf.n+...
                                (ai_m(sel(1))-1)*gridP.af.n*gridP.wf.n+...
                                (wi_new-1)*gridP.af.n+...
                                ai_new;
                            AJ_J(i,i) = AJ_J(i,i)-param.m*solutionS_s.el(ai_new,wi_new)^1*(ai_m_w(sel(1))*wi_m_w(sel(2)))/singles;
                            AJ_JJ(i,i_self)=AJ_JJ(i,i_self)+param.m*solutionS_s.el(ai_new,wi_new)^1*(ai_m_w(sel(1))*wi_m_w(sel(2)))/singles;
                        end                        
                    end
                    if always_marry||((solutionSP.JmlValue(ai, wi)<JJcValue)&&...
                       (solutionSP.JfcValue((ai_new-1)*param.Sdimmult+1, (wi_new-1)*param.Sdimmult+1)<JJcValue))
                        for k=1:4
                            sel='00';
                            binary = dec2bin(k-1);
                            sel(end-numel(binary)+1:end)=binary;
                            sel=reshape(str2num(char(num2cell(sel))),1,[]);
                            sel=sel+1;
                            i_self= (wi_m(sel(2))-1)*gridP.am.n*gridP.af.n*gridP.wf.n+...
                                (ai_m(sel(1))-1)*gridP.af.n*gridP.wf.n+...
                                (wi_new-1)*gridP.af.n+...
                                ai_new;
                            AJ_J(i,i) = AJ_J(i,i)-param.m*solutionS_s.ec(ai_new,wi_new)^1*(ai_m_w(sel(1))*wi_m_w(sel(2)))/singles;
                            AJ_JJc(i,i_self)=AJ_JJc(i,i_self)+param.m*solutionS_s.ec(ai_new,wi_new)^1*(ai_m_w(sel(1))*wi_m_w(sel(2)))/singles;
                        end   
                    end
                    if wi_new==1

                        UJValue=0;
                        UJcValue=0;
                        for k=1:4
                            sel='00';
                            binary = dec2bin(k-1);
                            sel(end-numel(binary)+1:end)=binary;
                            sel=reshape(str2num(char(num2cell(sel))),1,[]);
                            sel=sel+1;
                            UJValue=UJValue+...
                                solutionSP.UJlValue(ai_new, ai_m(sel(1)), wi_m(sel(2)))*...
                                (ai_m_w(sel(1))*wi_m_w(sel(2)));                            
                            UJcValue=UJcValue+...
                                solutionSP.UJcValue(ai_new, ai_m(sel(1)), wi_m(sel(2)))*...
                                (ai_m_w(sel(1))*wi_m_w(sel(2)));                            
                        end


                        if always_marry||((solutionSP.JmlValue(ai, wi)<UJValue)&&...
                           (solutionSP.UflValue((ai_new-1)*param.Sdimmult+1)<UJValue))
                            for k=1:4
                                sel='00';
                                binary = dec2bin(k-1);
                                sel(end-numel(binary)+1:end)=binary;
                                sel=reshape(str2num(char(num2cell(sel))),1,[]);
                                sel=sel+1;                  
                                i_self= (wi_m(sel(2))-1)*gridP.am.n*gridP.af.n+...
                                    (ai_m(sel(1))-1)*gridP.af.n+...
                                    ai_new;
                                AJ_J(i,i) = AJ_J(i,i)-param.m*solutionS_s.hl(ai_new)^1*(ai_m_w(sel(1))*wi_m_w(sel(2)))/singles;
                                AJ_UJ(i,i_self)=AJ_UJ(i,i_self)+param.m*solutionS_s.hl(ai_new)^1*(ai_m_w(sel(1))*wi_m_w(sel(2)))/singles;
                            end  

                        end
                        if always_marry||((solutionSP.JmlValue(ai, wi)<UJcValue)&&...
                           (solutionSP.UfcValue((ai_new-1)*param.Sdimmult+1)<UJcValue))
                            for k=1:4
                                sel='00';
                                binary = dec2bin(k-1);
                                sel(end-numel(binary)+1:end)=binary;
                                sel=reshape(str2num(char(num2cell(sel))),1,[]);
                                sel=sel+1;                  
                                i_self= (wi_m(sel(2))-1)*gridP.am.n*gridP.af.n+...
                                    (ai_m(sel(1))-1)*gridP.af.n+...
                                    ai_new;
                                AJ_J(i,i) = AJ_J(i,i)-param.m*solutionS_s.hc(ai_new)^1*(ai_m_w(sel(1))*wi_m_w(sel(2)))/singles;
                                AJ_UJc(i,i_self)=AJ_UJc(i,i_self)+param.m*solutionS_s.hc(ai_new)^1*(ai_m_w(sel(1))*wi_m_w(sel(2)))/singles;
                            end      
                        end

                        OJValue=0;
                        OJcValue=0;
                        for k=1:4
                            sel='00';
                            binary = dec2bin(k-1);
                            sel(end-numel(binary)+1:end)=binary;
                            sel=reshape(str2num(char(num2cell(sel))),1,[]);
                            sel=sel+1;
                            OJValue=OJValue+...
                                solutionSP.OJlValue(ai_new, ai_m(sel(1)), wi_m(sel(2)))*...
                                (ai_m_w(sel(1))*wi_m_w(sel(2)));                            
                            OJcValue=OJcValue+...
                                solutionSP.OJcValue(ai_new, ai_m(sel(1)), wi_m(sel(2)))*...
                                (ai_m_w(sel(1))*wi_m_w(sel(2)));                            
                        end
                        if always_marry||((solutionSP.JmlValue(ai, wi)<OJValue)&&...
                           (solutionSP.OflValue((ai_new-1)*param.Sdimmult+1)<OJValue))
                            for k=1:4
                                sel='00';
                                binary = dec2bin(k-1);
                                sel(end-numel(binary)+1:end)=binary;
                                sel=reshape(str2num(char(num2cell(sel))),1,[]);
                                sel=sel+1;                  
                                i_self= (wi_m(sel(2))-1)*gridP.am.n*gridP.af.n+...
                                    (ai_m(sel(1))-1)*gridP.af.n+...
                                    ai_new;
                                AJ_J(i,i) = AJ_J(i,i)-param.m*solutionS_s.ol(ai_new)^1*(ai_m_w(sel(1))*wi_m_w(sel(2)))/singles;
                                AJ_OJ(i,i_self)=AJ_OJ(i,i_self)+param.m*solutionS_s.ol(ai_new)^1*(ai_m_w(sel(1))*wi_m_w(sel(2)))/singles;
                            end  

                        end
                        if always_marry||((solutionSP.JmlValue(ai, wi)<OJcValue)&&...
                           (solutionSP.OfcValue((ai_new-1)*param.Sdimmult+1)<OJcValue))
                            for k=1:4
                                sel='00';
                                binary = dec2bin(k-1);
                                sel(end-numel(binary)+1:end)=binary;
                                sel=reshape(str2num(char(num2cell(sel))),1,[]);
                                sel=sel+1;                  
                                i_self= (wi_m(sel(2))-1)*gridP.am.n*gridP.af.n+...
                                    (ai_m(sel(1))-1)*gridP.af.n+...
                                    ai_new;
                                AJ_J(i,i) = AJ_J(i,i)-param.m*solutionS_s.oc(ai_new)^1*(ai_m_w(sel(1))*wi_m_w(sel(2)))/singles;
                                AJ_OJc(i,i_self)=AJ_OJc(i,i_self)+param.m*solutionS_s.oc(ai_new)^1*(ai_m_w(sel(1))*wi_m_w(sel(2)))/singles;
                            end      
                        end
                    end                
                end
            end
        end
    end   
    AJ_J = sparse(AJ_J)+pg.CJ;
    AJ_U = pg.CJ_U;  
    AJ_O = sparse(AJ_O);
    AJ_JJ = sparse(AJ_JJ);
    AJ_UJ = sparse(AJ_UJ);
    AJ_OJ = sparse(AJ_OJ);
    AJ_JJc = sparse(AJ_JJc);
    AJ_UJc = sparse(AJ_UJc);
    AJ_OJc = sparse(AJ_OJc);
end
function [AU_J, AU_U, AU_O, AU_UJ, AU_UU, AU_UO, AU_UJc, AU_UUc, AU_UOc] = generateAUl_f(solutionSP,solutionS_s, pg, gridS)
    AU_J = (zeros(gridS.a.n,pg.NS));
    AU_U = (zeros(gridS.a.n,gridS.a.n));
    AU_O = zeros(size(AU_U));
    AU_UU = (zeros(gridS.a.n,gridP.af.n*gridP.am.n));
    AU_UO = (zeros(gridS.a.n,gridP.af.n*gridP.am.n));
    AU_UJ = (zeros(gridS.a.n,param.NP/gridP.wf.n));
    AU_UUc = (zeros(gridS.a.n,gridP.af.n*gridP.am.n));
    AU_UOc = (zeros(gridS.a.n,gridP.af.n*gridP.am.n));
    AU_UJc = (zeros(gridS.a.n,param.NP/gridP.wf.n));
    singles=1;%sum(solutionSP.h,'all')+sum(solutionSP.e,'all');
    for ai = 1:gridS.a.n
        %find weights for linear interpolation
        ai_m=fVai_m(ai,:);
        ai_m_w=fVai_m_w(ai,:);
        if solutionSP.UflValue(ai) < solutionSP.OflValue(ai)
            AU_U(ai,ai) = AU_U(ai,ai) - pg.lambda_o;
            AU_O(ai,ai) = AU_O(ai,ai) + pg.lambda_o;
        elseif tremble~=0
            t = tremble/exp( solutionSP.UflValue(ai) - solutionSP.OflValue(ai) );
            AU_U(ai,ai) = AU_U(ai,ai) - t;
            AU_O(ai,ai) = AU_O(ai,ai) + t;
        end

        for wi_new = 1:gridS.w.n                
            i_new = (wi_new-1)*gridS.a.n + ai;
            if true&&(solutionSP.UflValue(ai)<solutionSP.JflValue(ai, wi_new))%solutionSP.UflValue(ai)<=(solutionSP.JflValue(ai, wi_new)-param.crit)
                AU_U(ai,ai) = AU_U(ai,ai)-pg.lambda_u*param.f_long(wi_new);
                AU_J(ai,i_new) = AU_J(ai,i_new)+pg.lambda_u*param.f_long(wi_new);
            elseif tremble~=0
                t = tremble/exp(solutionSP.UflValue(ai)-solutionSP.JflValue(ai, wi_new));
                AU_U(ai,ai) = AU_U(ai,ai)-t*param.f_long(wi_new);
                AU_J(ai,i_new) = AU_J(ai,i_new)+t*param.f_long(wi_new);
            end
        end
        for wi_new = 1:gridP.wm.n     
            %check marrige
            for ai_new = 1:gridP.am.n                

                UJValue=0;                
                UJcValue=0;                
                for k=1:2
                    sel='0';
                    binary = dec2bin(k-1);
                    sel(end-numel(binary)+1:end)=binary;
                    sel=reshape(str2num(char(num2cell(sel))),1,[]);
                    sel=sel+1;
                    UJValue=UJValue+...
                        solutionSP.UJlValue(ai_m(sel(1)), ai_new, wi_new)*...
                        (ai_m_w(sel(1)));                    
                    UJcValue=UJcValue+...
                        solutionSP.UJcValue(ai_m(sel(1)), ai_new, wi_new)*...
                        (ai_m_w(sel(1)));                    
                end

                if always_marry||((solutionSP.UflValue(ai)<UJValue)&&...
                   (solutionSP.JmlValue((ai_new-1)*param.Sdimmult+1, (wi_new-1)*param.Sdimmult+1)<UJValue))
                    for k=1:2
                        sel='0';
                        binary = dec2bin(k-1);
                        sel(end-numel(binary)+1:end)=binary;
                        sel=reshape(str2num(char(num2cell(sel))),1,[]);
                        sel=sel+1;
                        i_self = (wi_new-1)*gridP.af.n*gridP.am.n+(ai_new-1)*gridP.af.n+ai_m(sel(1));      
                        AU_U(ai,ai) = AU_U(ai,ai)-param.m*solutionS_s.el(ai_new,wi_new)^1*(ai_m_w(sel(1)))/singles;
                        AU_UJ(ai,i_self)=AU_UJ(ai,i_self)+param.m*solutionS_s.el(ai_new,wi_new)^1*(ai_m_w(sel(1)))/singles;
                    end                        
                end
                if always_marry||((solutionSP.UflValue(ai)<UJcValue)&&...
                   (solutionSP.JmcValue((ai_new-1)*param.Sdimmult+1, (wi_new-1)*param.Sdimmult+1)<UJcValue))
                    for k=1:2
                        sel='0';
                        binary = dec2bin(k-1);
                        sel(end-numel(binary)+1:end)=binary;
                        sel=reshape(str2num(char(num2cell(sel))),1,[]);
                        sel=sel+1;
                        i_self = (wi_new-1)*gridP.af.n*gridP.am.n+(ai_new-1)*gridP.af.n+ai_m(sel(1));      
                        AU_U(ai,ai) = AU_U(ai,ai)-param.m*solutionS_s.ec(ai_new,wi_new)^1*(ai_m_w(sel(1)))/singles;
                        AU_UJc(ai,i_self)=AU_UJc(ai,i_self)+param.m*solutionS_s.ec(ai_new,wi_new)^1*(ai_m_w(sel(1)))/singles;
                    end   
                end
                if wi_new==1

                    UUValue=0;
                    UUcValue=0;
                    for k=1:2
                        sel='0';
                        binary = dec2bin(k-1);
                        sel(end-numel(binary)+1:end)=binary;
                        sel=reshape(str2num(char(num2cell(sel))),1,[]);
                        sel=sel+1;
                        UUValue=UUValue+...
                            solutionSP.UUlValue(ai_m(sel(1)), ai_new)*...
                            (ai_m_w(sel(1)));                        
                        UUcValue=UUcValue+...
                            solutionSP.UUcValue(ai_m(sel(1)), ai_new)*...
                            (ai_m_w(sel(1)));                        
                    end


                    if always_marry||((solutionSP.UflValue(ai)<UUValue)&&...
                       (solutionSP.UmlValue((ai_new-1)*param.Sdimmult+1)<UUValue))
                        for k=1:2
                            sel='0';
                            binary = dec2bin(k-1);
                            sel(end-numel(binary)+1:end)=binary;
                            sel=reshape(str2num(char(num2cell(sel))),1,[]);
                            sel=sel+1;                  
                            i_self= (ai_new-1)*gridP.af.n+ai_m(sel(1));
                            AU_U(ai,ai) = AU_U(ai,ai)-param.m*solutionS_s.hl(ai_new)^1*(ai_m_w(sel(1)))/singles;
                            AU_UU(ai,i_self)=AU_UU(ai,i_self)+param.m*solutionS_s.hl(ai_new)^1*(ai_m_w(sel(1)))/singles;
                        end                          
                    end
                    if always_marry||((solutionSP.UflValue(ai)<UUcValue)&&...
                       (solutionSP.UmcValue((ai_new-1)*param.Sdimmult+1)<UUcValue))
                        for k=1:2
                            sel='0';
                            binary = dec2bin(k-1);
                            sel(end-numel(binary)+1:end)=binary;
                            sel=reshape(str2num(char(num2cell(sel))),1,[]);
                            sel=sel+1;                  
                            i_self= (ai_new-1)*gridP.af.n+ai_m(sel(1));
                            AU_U(ai,ai) = AU_U(ai,ai)-param.m*solutionS_s.hc(ai_new)^1*(ai_m_w(sel(1)))/singles;
                            AU_UUc(ai,i_self)=AU_UUc(ai,i_self)+param.m*solutionS_s.hc(ai_new)^1*(ai_m_w(sel(1)))/singles;
                        end     
                    end

                    UOValue=0;
                    UOcValue=0;
                    for k=1:2
                        sel='0';
                        binary = dec2bin(k-1);
                        sel(end-numel(binary)+1:end)=binary;
                        sel=reshape(str2num(char(num2cell(sel))),1,[]);
                        sel=sel+1;
                        UOValue=UOValue+...
                            solutionSP.UOlValue(ai_m(sel(1)), ai_new)*...
                            (ai_m_w(sel(1)));                        
                        UOcValue=UOcValue+...
                            solutionSP.UOcValue(ai_m(sel(1)), ai_new)*...
                            (ai_m_w(sel(1)));                        
                    end
                    if always_marry||((solutionSP.UflValue(ai)<UOValue)&&...
                       (solutionSP.OmlValue((ai_new-1)*param.Sdimmult+1)<UOValue))
                        for k=1:2
                            sel='0';
                            binary = dec2bin(k-1);
                            sel(end-numel(binary)+1:end)=binary;
                            sel=reshape(str2num(char(num2cell(sel))),1,[]);
                            sel=sel+1;                  
                            i_self= (ai_new-1)*gridP.af.n+ai_m(sel(1));
                            AU_U(ai,ai) = AU_U(ai,ai)-param.m*solutionS_s.ol(ai_new)^1*(ai_m_w(sel(1)))/singles;
                            AU_UO(ai,i_self)=AU_UO(ai,i_self)+param.m*solutionS_s.ol(ai_new)^1*(ai_m_w(sel(1)))/singles;
                        end                          
                    end
                    if always_marry||((solutionSP.UflValue(ai)<UOcValue)&&...
                       (solutionSP.OmcValue((ai_new-1)*param.Sdimmult+1)<UOcValue))
                        for k=1:2
                            sel='0';
                            binary = dec2bin(k-1);
                            sel(end-numel(binary)+1:end)=binary;
                            sel=reshape(str2num(char(num2cell(sel))),1,[]);
                            sel=sel+1;                  
                            i_self= (ai_new-1)*gridP.af.n+ai_m(sel(1));
                            AU_U(ai,ai) = AU_U(ai,ai)-param.m*solutionS_s.oc(ai_new)^1*(ai_m_w(sel(1)))/singles;
                            AU_UOc(ai,i_self)=AU_UOc(ai,i_self)+param.m*solutionS_s.oc(ai_new)^1*(ai_m_w(sel(1)))/singles;
                        end     
                    end
                end

            end
        end
    end 
    AU_J = sparse(AU_J);
    AU_U = sparse(AU_U) + pg.CU;
    AU_O = sparse(AU_O);
    AU_UU = sparse(AU_UU);
    AU_UO = sparse(AU_UO);
    AU_UJ = sparse(AU_UJ);
    AU_UUc = sparse(AU_UUc);
    AU_UOc = sparse(AU_UOc);
    AU_UJc = sparse(AU_UJc);
end
function [AU_J, AU_U, AU_O, AU_JU, AU_UU, AU_OU, AU_JUc, AU_UUc, AU_OUc] = generateAUl_m(solutionSP,solutionS_s, pg, gridS)
    AU_J = (zeros(gridS.a.n,pg.NS));
    AU_U = (zeros(gridS.a.n,gridS.a.n));
    AU_O = zeros(size(AU_U));
    AU_UU = (zeros(gridS.a.n,gridP.af.n*gridP.am.n));
    AU_OU = (zeros(gridS.a.n,gridP.af.n*gridP.am.n));
    AU_JU = (zeros(gridS.a.n,param.NP/gridP.wm.n));
    AU_UUc = (zeros(gridS.a.n,gridP.af.n*gridP.am.n));
    AU_OUc = (zeros(gridS.a.n,gridP.af.n*gridP.am.n));
    AU_JUc = (zeros(gridS.a.n,param.NP/gridP.wm.n));
    singles=1;%sum(solutionSP.h,'all')+sum(solutionSP.e,'all');
    for ai = 1:gridS.a.n
        %find weights for linear interpolation
        ai_m=mVai_m(ai,:);
        ai_m_w=mVai_m_w(ai,:);
        if solutionSP.UmlValue(ai) < solutionSP.OmlValue(ai)
            AU_U(ai,ai) = AU_U(ai,ai) - pg.lambda_o;
            AU_O(ai,ai) = AU_O(ai,ai) + pg.lambda_o;
        elseif tremble~=0
            t = tremble/exp( solutionSP.UmlValue(ai) - solutionSP.OmlValue(ai) );
            AU_U(ai,ai) = AU_U(ai,ai) - t;
            AU_O(ai,ai) = AU_O(ai,ai) + t;
        end
        for wi_new = 1:gridS.w.n                
            i_new = (wi_new-1)*gridS.a.n + ai;
            if true&&(solutionSP.UmlValue(ai)<solutionSP.JmlValue(ai, wi_new))%solutionSP.UmlValue(ai)<=(solutionSP.JmlValue(ai, wi_new)-param.crit)
                AU_U(ai,ai) = AU_U(ai,ai)-pg.lambda_u*param.f_long(wi_new);
                AU_J(ai,i_new) = AU_J(ai,i_new)+pg.lambda_u*param.f_long(wi_new);
            elseif tremble~=0
                t = tremble/exp(solutionSP.UmlValue(ai)-solutionSP.JmlValue(ai, wi_new));
                AU_U(ai,ai) = AU_U(ai,ai)-t*param.f_long(wi_new);
                AU_J(ai,i_new) = AU_J(ai,i_new)+t*param.f_long(wi_new);
            end
        end
        for wi_new = 1:gridP.wf.n    % was wm 
            %check marrige
            for ai_new = 1:gridP.am.n                

                JUValue=0;                
                JUcValue=0;                
                for k=1:2
                    sel='0';
                    binary = dec2bin(k-1);
                    sel(end-numel(binary)+1:end)=binary;
                    sel=reshape(str2num(char(num2cell(sel))),1,[]);
                    sel=sel+1;
                    JUValue=JUValue+...
                        solutionSP.JUlValue(ai_new, wi_new, ai_m(sel(1)))*...
                        (ai_m_w(sel(1)));                    
                    JUcValue=JUcValue+...
                        solutionSP.JUcValue(ai_new, wi_new, ai_m(sel(1)))*...
                        (ai_m_w(sel(1)));                    
                end

                if always_marry||((solutionSP.UmlValue(ai)<JUValue)&&...
                   (solutionSP.JflValue((ai_new-1)*param.Sdimmult+1, (wi_new-1)*param.Sdimmult+1)<JUValue))
                    for k=1:2
                        sel='0';
                        binary = dec2bin(k-1);
                        sel(end-numel(binary)+1:end)=binary;
                        sel=reshape(str2num(char(num2cell(sel))),1,[]);
                        sel=sel+1;             
                        i_self= (ai_m(sel(1))-1)*gridP.af.n*gridP.wf.n+...
                                (wi_new-1)*gridP.af.n+...
                                ai_new;
                        AU_U(ai,ai) = AU_U(ai,ai)-param.m*solutionS_s.el(ai_new,wi_new)^1*(ai_m_w(sel(1)))/singles;
                        AU_JU(ai,i_self)=AU_JU(ai,i_self)+param.m*solutionS_s.el(ai_new,wi_new)^1*(ai_m_w(sel(1)))/singles;
                    end                        
                end
                if always_marry||((solutionSP.UmlValue(ai)<JUcValue)&&...
                   (solutionSP.JfcValue((ai_new-1)*param.Sdimmult+1, (wi_new-1)*param.Sdimmult+1)<JUcValue))
                    for k=1:2
                        sel='0';
                        binary = dec2bin(k-1);
                        sel(end-numel(binary)+1:end)=binary;
                        sel=reshape(str2num(char(num2cell(sel))),1,[]);
                        sel=sel+1;
                        i_self= (ai_m(sel(1))-1)*gridP.af.n*gridP.wf.n+...
                                (wi_new-1)*gridP.af.n+...
                                ai_new;      
                        AU_U(ai,ai) = AU_U(ai,ai)-param.m*solutionS_s.ec(ai_new,wi_new)^1*(ai_m_w(sel(1)))/singles;
                        AU_JUc(ai,i_self)=AU_JUc(ai,i_self)+param.m*solutionS_s.ec(ai_new,wi_new)^1*(ai_m_w(sel(1)))/singles;
                    end   
                end
                if wi_new==1

                    UUValue=0;
                    UUcValue=0;
                    for k=1:2
                        sel='0';
                        binary = dec2bin(k-1);
                        sel(end-numel(binary)+1:end)=binary;
                        sel=reshape(str2num(char(num2cell(sel))),1,[]);
                        sel=sel+1;
                        UUValue=UUValue+...
                            solutionSP.UUlValue(ai_m(sel(1)), ai_new)*...
                            (ai_m_w(sel(1)));                        
                        UUcValue=UUcValue+...
                            solutionSP.UUcValue(ai_m(sel(1)), ai_new)*...
                            (ai_m_w(sel(1)));                        
                    end


                    if always_marry||((solutionSP.UmlValue(ai)<UUValue)&&...
                       (solutionSP.UflValue((ai_new-1)*param.Sdimmult+1)<UUValue))
                        for k=1:2
                            sel='0';
                            binary = dec2bin(k-1);
                            sel(end-numel(binary)+1:end)=binary;
                            sel=reshape(str2num(char(num2cell(sel))),1,[]);
                            sel=sel+1;                  
                            i_self= (ai_m(sel(1))-1)*gridP.af.n+ai_new;
                            AU_U(ai,ai) = AU_U(ai,ai)-param.m*solutionS_s.hl(ai_new)^1*(ai_m_w(sel(1)))/singles;
                            AU_UU(ai,i_self)=AU_UU(ai,i_self)+param.m*solutionS_s.hl(ai_new)^1*(ai_m_w(sel(1)))/singles;
                        end                          
                    end
                    if always_marry||((solutionSP.UmlValue(ai)<UUcValue)&&...
                       (solutionSP.UfcValue((ai_new-1)*param.Sdimmult+1)<UUcValue))
                        for k=1:2
                            sel='0';
                            binary = dec2bin(k-1);
                            sel(end-numel(binary)+1:end)=binary;
                            sel=reshape(str2num(char(num2cell(sel))),1,[]);
                            sel=sel+1;                  
                            i_self= (ai_m(sel(1))-1)*gridP.af.n+ai_new;
                            AU_U(ai,ai) = AU_U(ai,ai)-param.m*solutionS_s.hc(ai_new)^1*(ai_m_w(sel(1)))/singles;
                            AU_UUc(ai,i_self)=AU_UUc(ai,i_self)+param.m*solutionS_s.hc(ai_new)^1*(ai_m_w(sel(1)))/singles;
                        end     
                    end

                    OUValue=0;
                    OUcValue=0;
                    for k=1:2
                        sel='0';
                        binary = dec2bin(k-1);
                        sel(end-numel(binary)+1:end)=binary;
                        sel=reshape(str2num(char(num2cell(sel))),1,[]);
                        sel=sel+1;
                        OUValue=OUValue+...
                            solutionSP.OUlValue(ai_m(sel(1)), ai_new)*...
                            (ai_m_w(sel(1)));                        
                        OUcValue=OUcValue+...
                            solutionSP.OUcValue(ai_m(sel(1)), ai_new)*...
                            (ai_m_w(sel(1)));                        
                    end
                    if always_marry||((solutionSP.UmlValue(ai)<OUValue)&&...
                       (solutionSP.OflValue((ai_new-1)*param.Sdimmult+1)<OUValue))
                        for k=1:2
                            sel='0';
                            binary = dec2bin(k-1);
                            sel(end-numel(binary)+1:end)=binary;
                            sel=reshape(str2num(char(num2cell(sel))),1,[]);
                            sel=sel+1;                  
                            i_self= (ai_m(sel(1))-1)*gridP.af.n+ai_new;
                            AU_U(ai,ai) = AU_U(ai,ai)-param.m*solutionS_s.ol(ai_new)^1*(ai_m_w(sel(1)))/singles;
                            AU_OU(ai,i_self)=AU_OU(ai,i_self)+param.m*solutionS_s.ol(ai_new)^1*(ai_m_w(sel(1)))/singles;
                        end                          
                    end
                    if always_marry||((solutionSP.UmlValue(ai)<OUcValue)&&...
                       (solutionSP.OfcValue((ai_new-1)*param.Sdimmult+1)<OUcValue))
                        for k=1:2
                            sel='0';
                            binary = dec2bin(k-1);
                            sel(end-numel(binary)+1:end)=binary;
                            sel=reshape(str2num(char(num2cell(sel))),1,[]);
                            sel=sel+1;                  
                            i_self= (ai_m(sel(1))-1)*gridP.af.n+ai_new;
                            AU_U(ai,ai) = AU_U(ai,ai)-param.m*solutionS_s.oc(ai_new)^1*(ai_m_w(sel(1)))/singles;
                            AU_OUc(ai,i_self)=AU_OUc(ai,i_self)+param.m*solutionS_s.oc(ai_new)^1*(ai_m_w(sel(1)))/singles;
                        end     
                    end
                end

            end
        end
    end 
    AU_J = sparse(AU_J);
    AU_U = sparse(AU_U) + pg.CU;
    AU_O = sparse(AU_O);
    AU_UU = sparse(AU_UU);
    AU_OU = sparse(AU_OU);
    AU_JU = sparse(AU_JU);
    AU_UUc = sparse(AU_UUc);
    AU_OUc = sparse(AU_OUc);
    AU_JUc = sparse(AU_JUc);
end
function [AO_O, AO_J, AO_OJ, AO_OU, AO_OJc, AO_OUc] = generateAOl_f(solutionSP,solutionS_s, pg, gridS)    
    AO_O  = (zeros(gridS.a.n,gridS.a.n));
    AO_J  = (zeros(gridS.a.n,pg.NS));              % NEW: O -> J block (dest is single-J grid)
    AO_OU = (zeros(gridS.a.n,gridP.af.n*gridP.am.n));
    AO_OJ = (zeros(gridS.a.n,param.NP/gridP.wf.n));
    AO_OUc= (zeros(gridS.a.n,gridP.af.n*gridP.am.n));
    AO_OJc= (zeros(gridS.a.n,param.NP/gridP.wf.n));
    singles=1;

    for ai = 1:gridS.a.n
        ai_m   = fVai_m(ai,:);
        ai_m_w = fVai_m_w(ai,:);

        % ---- O -> J over wage grid (replaces O -> U) ----
        for wi_new = 1:gridS.w.n
            i_new = (wi_new-1)*gridS.a.n + ai;   % single-J index (w,a)
            if solutionSP.OflValue(ai) < solutionSP.JflValue(ai, wi_new)
                AO_O(ai,ai) = AO_O(ai,ai) - param.pf.lambda_oE * param.f_long(wi_new);
                AO_J(ai,i_new) = AO_J(ai,i_new) + param.pf.lambda_oE * param.f_long(wi_new);
            elseif tremble~=0
                t = tremble/exp( solutionSP.OflValue(ai) - solutionSP.JflValue(ai, wi_new) );
                AO_O(ai,ai) = AO_O(ai,ai) - t * param.f_long(wi_new);
                AO_J(ai,i_new) = AO_J(ai,i_new) + t * param.f_long(wi_new);
            end
        end
        for wi_new = 1:gridP.wm.n   
            %check marrige
            for ai_new = 1:gridP.am.n                

                OJValue=0;                
                OJcValue=0;                
                for k=1:2
                    sel='0';
                    binary = dec2bin(k-1);
                    sel(end-numel(binary)+1:end)=binary;
                    sel=reshape(str2num(char(num2cell(sel))),1,[]);
                    sel=sel+1;
                    OJValue=OJValue+...
                        solutionSP.OJlValue(ai_m(sel(1)), ai_new, wi_new)*...
                        (ai_m_w(sel(1)));                    
                    OJcValue=OJcValue+...
                        solutionSP.OJcValue(ai_m(sel(1)), ai_new, wi_new)*...
                        (ai_m_w(sel(1)));                    
                end

                if always_marry||((solutionSP.OflValue(ai)<OJValue)&&...
                   (solutionSP.JmlValue((ai_new-1)*param.Sdimmult+1, (wi_new-1)*param.Sdimmult+1)<OJValue))
                    for k=1:2
                        sel='0';
                        binary = dec2bin(k-1);
                        sel(end-numel(binary)+1:end)=binary;
                        sel=reshape(str2num(char(num2cell(sel))),1,[]);
                        sel=sel+1;
                        i_self = (wi_new-1)*gridP.af.n*gridP.am.n+(ai_new-1)*gridP.af.n+ai_m(sel(1));      
                        AO_O(ai,ai) = AO_O(ai,ai)-param.m*solutionS_s.el(ai_new,wi_new)^1*(ai_m_w(sel(1)))/singles;
                        AO_OJ(ai,i_self)=AO_OJ(ai,i_self)+param.m*solutionS_s.el(ai_new,wi_new)^1*(ai_m_w(sel(1)))/singles;
                    end                        
                end
                if always_marry||((solutionSP.OflValue(ai)<OJcValue)&&...
                   (solutionSP.JmcValue((ai_new-1)*param.Sdimmult+1, (wi_new-1)*param.Sdimmult+1)<OJcValue))
                    for k=1:2
                        sel='0';
                        binary = dec2bin(k-1);
                        sel(end-numel(binary)+1:end)=binary;
                        sel=reshape(str2num(char(num2cell(sel))),1,[]);
                        sel=sel+1;
                        i_self = (wi_new-1)*gridP.af.n*gridP.am.n+(ai_new-1)*gridP.af.n+ai_m(sel(1));      
                        AO_O(ai,ai) = AO_O(ai,ai)-param.m*solutionS_s.ec(ai_new,wi_new)^1*(ai_m_w(sel(1)))/singles;
                        AO_OJc(ai,i_self)=AO_OJc(ai,i_self)+param.m*solutionS_s.ec(ai_new,wi_new)^1*(ai_m_w(sel(1)))/singles;
                    end   
                end
                if wi_new==1

                    OUValue=0;
                    OUcValue=0;
                    for k=1:2
                        sel='0';
                        binary = dec2bin(k-1);
                        sel(end-numel(binary)+1:end)=binary;
                        sel=reshape(str2num(char(num2cell(sel))),1,[]);
                        sel=sel+1;
                        OUValue=OUValue+...
                            solutionSP.OUlValue(ai_m(sel(1)), ai_new)*...
                            (ai_m_w(sel(1)));                        
                        OUcValue=OUcValue+...
                            solutionSP.OUcValue(ai_m(sel(1)), ai_new)*...
                            (ai_m_w(sel(1)));                        
                    end


                    if always_marry||((solutionSP.OflValue(ai)<OUValue)&&...
                       (solutionSP.UmlValue((ai_new-1)*param.Sdimmult+1)<OUValue))
                        for k=1:2
                            sel='0';
                            binary = dec2bin(k-1);
                            sel(end-numel(binary)+1:end)=binary;
                            sel=reshape(str2num(char(num2cell(sel))),1,[]);
                            sel=sel+1;                  
                            i_self= (ai_new-1)*gridP.af.n+ai_m(sel(1));
                            AO_O(ai,ai) = AO_O(ai,ai)-param.m*solutionS_s.hl(ai_new)^1*(ai_m_w(sel(1)))/singles;
                            AO_OU(ai,i_self)=AO_OU(ai,i_self)+param.m*solutionS_s.hl(ai_new)^1*(ai_m_w(sel(1)))/singles;
                        end                          
                    end
                    if always_marry||((solutionSP.OflValue(ai)<OUcValue)&&...
                       (solutionSP.UmcValue((ai_new-1)*param.Sdimmult+1)<OUcValue))
                        for k=1:2
                            sel='0';
                            binary = dec2bin(k-1);
                            sel(end-numel(binary)+1:end)=binary;
                            sel=reshape(str2num(char(num2cell(sel))),1,[]);
                            sel=sel+1;                  
                            i_self= (ai_new-1)*gridP.af.n+ai_m(sel(1));
                            AO_O(ai,ai) = AO_O(ai,ai)-param.m*solutionS_s.hc(ai_new)^1*(ai_m_w(sel(1)))/singles;
                            AO_OUc(ai,i_self)=AO_OUc(ai,i_self)+param.m*solutionS_s.hc(ai_new)^1*(ai_m_w(sel(1)))/singles;
                        end     
                    end
                end

            end
        end
    end     
    AO_O = sparse(AO_O) + pg.CO;
    AO_J = sparse(AO_J);
    AO_OU = sparse(AO_OU);
    AO_OJ = sparse(AO_OJ);
    AO_OUc = sparse(AO_OUc);
    AO_OJc = sparse(AO_OJc);
end

function [AO_O, AO_J, AO_JO, AO_UO, AO_JOc, AO_UOc] = generateAOl_m(solutionSP,solutionS_s, pg, gridS)    
    AO_O  = (zeros(gridS.a.n,gridS.a.n));
    AO_J  = (zeros(gridS.a.n,pg.NS));              % NEW: O -> J block (dest is single-J grid)
    AO_UO = (zeros(gridS.a.n,gridP.af.n*gridP.am.n));
    AO_JO = (zeros(gridS.a.n,param.NP/gridP.wm.n));
    AO_UOc= (zeros(gridS.a.n,gridP.af.n*gridP.am.n));
    AO_JOc= (zeros(gridS.a.n,param.NP/gridP.wm.n));
    singles=1;

    for ai = 1:gridS.a.n
        ai_m   = mVai_m(ai,:);
        ai_m_w = mVai_m_w(ai,:);

        % ---- O -> J over wage grid (replaces O -> U) ----
        for wi_new = 1:gridS.w.n
            i_new = (wi_new-1)*gridS.a.n + ai;   % single-J index (w,a)
            if solutionSP.OmlValue(ai) < solutionSP.JmlValue(ai, wi_new)
                AO_O(ai,ai) = AO_O(ai,ai) - param.pm.lambda_oE * param.f_long(wi_new);
                AO_J(ai,i_new) = AO_J(ai,i_new) + param.pm.lambda_oE * param.f_long(wi_new);
            elseif tremble~=0
                t = tremble/exp( solutionSP.OmlValue(ai) - solutionSP.JmlValue(ai, wi_new) );
                AO_O(ai,ai) = AO_O(ai,ai) - t * param.f_long(wi_new);
                AO_J(ai,i_new) = AO_J(ai,i_new) + t * param.f_long(wi_new);
            end
        end

        for wi_new = 1:gridP.wf.n % was wm     
            %check marrige
            for ai_new = 1:gridP.am.n                

                JOValue=0;                
                JOcValue=0;                
                for k=1:2
                    sel='0';
                    binary = dec2bin(k-1);
                    sel(end-numel(binary)+1:end)=binary;
                    sel=reshape(str2num(char(num2cell(sel))),1,[]);
                    sel=sel+1;
                    JOValue=JOValue+...
                        solutionSP.JOlValue(ai_new, wi_new, ai_m(sel(1)))*...
                        (ai_m_w(sel(1)));                    
                    JOcValue=JOcValue+...
                        solutionSP.JOcValue(ai_new, wi_new, ai_m(sel(1)))*...
                        (ai_m_w(sel(1)));                    
                end

                if always_marry||((solutionSP.OmlValue(ai)<JOValue)&&...
                   (solutionSP.JflValue((ai_new-1)*param.Sdimmult+1, (wi_new-1)*param.Sdimmult+1)<JOValue))
                    for k=1:2
                        sel='0';
                        binary = dec2bin(k-1);
                        sel(end-numel(binary)+1:end)=binary;
                        sel=reshape(str2num(char(num2cell(sel))),1,[]);
                        sel=sel+1;             
                        i_self= (ai_m(sel(1))-1)*gridP.af.n*gridP.wf.n+...
                                (wi_new-1)*gridP.af.n+...
                                ai_new;
                        AO_O(ai,ai) = AO_O(ai,ai)-param.m*solutionS_s.el(ai_new,wi_new)^1*(ai_m_w(sel(1)))/singles;
                        AO_JO(ai,i_self)=AO_JO(ai,i_self)+param.m*solutionS_s.el(ai_new,wi_new)^1*(ai_m_w(sel(1)))/singles;
                    end                        
                end
                if always_marry||((solutionSP.OmlValue(ai)<JOcValue)&&...
                   (solutionSP.JfcValue((ai_new-1)*param.Sdimmult+1, (wi_new-1)*param.Sdimmult+1)<JOcValue))
                    for k=1:2
                        sel='0';
                        binary = dec2bin(k-1);
                        sel(end-numel(binary)+1:end)=binary;
                        sel=reshape(str2num(char(num2cell(sel))),1,[]);
                        sel=sel+1;
                        i_self= (ai_m(sel(1))-1)*gridP.af.n*gridP.wf.n+...
                                (wi_new-1)*gridP.af.n+...
                                ai_new;      
                        AO_O(ai,ai) = AO_O(ai,ai)-param.m*solutionS_s.ec(ai_new,wi_new)^1*(ai_m_w(sel(1)))/singles;
                        AO_JOc(ai,i_self)=AO_JOc(ai,i_self)+param.m*solutionS_s.ec(ai_new,wi_new)^1*(ai_m_w(sel(1)))/singles;
                    end   
                end
                if wi_new==1

                    UOValue=0;
                    UOcValue=0;
                    for k=1:2
                        sel='0';
                        binary = dec2bin(k-1);
                        sel(end-numel(binary)+1:end)=binary;
                        sel=reshape(str2num(char(num2cell(sel))),1,[]);
                        sel=sel+1;
                        UOValue=UOValue+...
                            solutionSP.UOlValue(ai_new,ai_m(sel(1)))*...
                            (ai_m_w(sel(1)));                        
                        UOcValue=UOcValue+...
                            solutionSP.UOcValue(ai_new,ai_m(sel(1)))*...
                            (ai_m_w(sel(1)));                        
                    end


                    if always_marry||((solutionSP.OmlValue(ai)<UOValue)&&...
                       (solutionSP.UflValue((ai_new-1)*param.Sdimmult+1)<UOValue))
                        for k=1:2
                            sel='0';
                            binary = dec2bin(k-1);
                            sel(end-numel(binary)+1:end)=binary;
                            sel=reshape(str2num(char(num2cell(sel))),1,[]);
                            sel=sel+1;                  
                            i_self= (ai_m(sel(1))-1)*gridP.af.n+ai_new;
                            AO_O(ai,ai) = AO_O(ai,ai)-param.m*solutionS_s.hl(ai_new)^1*(ai_m_w(sel(1)))/singles;
                            AO_UO(ai,i_self)=AO_UO(ai,i_self)+param.m*solutionS_s.hl(ai_new)^1*(ai_m_w(sel(1)))/singles;
                        end                          
                    end
                    if always_marry||((solutionSP.OmlValue(ai)<UOcValue)&&...
                       (solutionSP.UfcValue((ai_new-1)*param.Sdimmult+1)<UOcValue))
                        for k=1:2
                            sel='0';
                            binary = dec2bin(k-1);
                            sel(end-numel(binary)+1:end)=binary;
                            sel=reshape(str2num(char(num2cell(sel))),1,[]);
                            sel=sel+1;                  
                            i_self= (ai_m(sel(1))-1)*gridP.af.n+ai_new;
                            AO_O(ai,ai) = AO_O(ai,ai)-param.m*solutionS_s.hc(ai_new)^1*(ai_m_w(sel(1)))/singles;
                            AO_UOc(ai,i_self)=AO_UOc(ai,i_self)+param.m*solutionS_s.hc(ai_new)^1*(ai_m_w(sel(1)))/singles;
                        end     
                    end
                end

            end
        end
    end     
    AO_O = sparse(AO_O) + pg.CO;
    AO_J = sparse(AO_J);
    AO_UO = sparse(AO_UO);
    AO_JO = sparse(AO_JO);
    AO_UOc = sparse(AO_UOc);
    AO_JOc = sparse(AO_JOc);
end

function [AJ_J, AJ_U, AJ_O, AJ_JJ, AJ_JU, AJ_JO] = generateAJc_f(solutionSP,solutionS_s, pg, gridS)
    AJ_J = (zeros(pg.NS,pg.NS));
    AJ_O = zeros(size(pg.CJ_U));
    AJ_JJ = (zeros(pg.NS,param.NP));
    AJ_JU = (zeros(pg.NS,param.NP/gridP.wm.n));
    AJ_JO = (zeros(pg.NS,param.NP/gridP.wm.n));
    singles=1;%sum(solutionSP.h,'all')+sum(solutionSP.e,'all');

    for ai = 1:gridS.a.n
        %find weights for linear interpolation
        ai_m=fVai_m(ai,:);
        ai_m_w=fVai_m_w(ai,:);

        for wi = 1:gridS.w.n
            %find weights for linear interpolation
            wi_m=fVwi_m(wi,:);
            wi_m_w=fVwi_m_w(wi,:);

            i = (wi-1)*gridS.a.n + ai;
            if solutionSP.JfcValue(ai, wi) < solutionSP.OfcValue(ai)
                AJ_J(i,i)   = AJ_J(i,i)   - pg.lambda_o;
                AJ_O(i,ai)  = AJ_O(i,ai)  + pg.lambda_o;
            elseif tremble~=0
                t = tremble/exp( solutionSP.JfcValue(ai, wi) - solutionSP.OfcValue(ai) );
                AJ_J(i,i)   = AJ_J(i,i)   - t;
                AJ_O(i,ai)  = AJ_O(i,ai)  + t;
            end
            % for wi_new = 1:gridS.w.n
            %     i_new = (wi_new-1)*gridS.a.n + ai;
            %     if solutionSP.JfcValue(ai, wi)<solutionSP.JfcValue(ai, wi_new)
            %         AJ_J(i,i) = AJ_J(i,i)-pg.lambda*param.f_long(wi_new);
            %         AJ_J(i,i_new) = AJ_J(i,i_new)+pg.lambda*param.f_long(wi_new);
            %     end
            % end
            for wi_new = 1:gridP.wm.n
                %====
                %check marrige
                for ai_new = 1:gridP.am.n                                   
                    JJcValue=0;                    
                    for k=1:4
                        sel='00';
                        binary = dec2bin(k-1);
                        sel(end-numel(binary)+1:end)=binary;
                        sel=reshape(str2num(char(num2cell(sel))),1,[]);
                        sel=sel+1;                                              
                        JJcValue=JJcValue+...
                            solutionSP.JJcValue(ai_m(sel(1)), wi_m(sel(2)), ai_new, wi_new)*...
                            (ai_m_w(sel(1))*wi_m_w(sel(2)));                        
                    end


                    if always_marry||((solutionSP.JfcValue(ai, wi)<JJcValue)&&...
                       (solutionSP.JmlValue((ai_new-1)*param.Sdimmult+1, (wi_new-1)*param.Sdimmult+1)<JJcValue))
                        for k=1:4
                            sel='00';
                            binary = dec2bin(k-1);
                            sel(end-numel(binary)+1:end)=binary;
                            sel=reshape(str2num(char(num2cell(sel))),1,[]);
                            sel=sel+1;
                            i_self= (wi_new-1)*gridP.am.n*gridP.af.n*gridP.wf.n+...
                                (ai_new-1)*gridP.af.n*gridP.wf.n+(wi_m(sel(2))-1)*gridP.af.n+ai_m(sel(1));
                            AJ_J(i,i) = AJ_J(i,i)-param.m*solutionS_s.el(ai_new,wi_new)^1*(ai_m_w(sel(1))*wi_m_w(sel(2)))/singles;
                            AJ_JJ(i,i_self)=AJ_JJ(i,i_self)+param.m*solutionS_s.el(ai_new,wi_new)^1*(ai_m_w(sel(1))*wi_m_w(sel(2)))/singles;
                        end                        
                    end
                    if always_marry||((solutionSP.JfcValue(ai, wi)<JJcValue)&&...
                       (solutionSP.JmcValue((ai_new-1)*param.Sdimmult+1, (wi_new-1)*param.Sdimmult+1)<JJcValue))
                        for k=1:4
                            sel='00';
                            binary = dec2bin(k-1);
                            sel(end-numel(binary)+1:end)=binary;
                            sel=reshape(str2num(char(num2cell(sel))),1,[]);
                            sel=sel+1;
                            i_self= (wi_new-1)*gridP.am.n*gridP.af.n*gridP.wf.n+...
                                (ai_new-1)*gridP.af.n*gridP.wf.n+(wi_m(sel(2))-1)*gridP.af.n+ai_m(sel(1));
                            AJ_J(i,i) = AJ_J(i,i)-param.m*solutionS_s.ec(ai_new,wi_new)^1*(ai_m_w(sel(1))*wi_m_w(sel(2)))/singles;
                            AJ_JJ(i,i_self)=AJ_JJ(i,i_self)+param.m*solutionS_s.ec(ai_new,wi_new)^1*(ai_m_w(sel(1))*wi_m_w(sel(2)))/singles;
                        end   
                    end
                    if wi_new==1

                        JUcValue=0;
                        for k=1:4
                            sel='00';
                            binary = dec2bin(k-1);
                            sel(end-numel(binary)+1:end)=binary;
                            sel=reshape(str2num(char(num2cell(sel))),1,[]);
                            sel=sel+1;                         
                            JUcValue=JUcValue+...
                                solutionSP.JUcValue(ai_m(sel(1)), wi_m(sel(2)), ai_new)*...
                                (ai_m_w(sel(1))*wi_m_w(sel(2)));                            
                        end


                        if always_marry||((solutionSP.JfcValue(ai, wi)<JUcValue)&&...
                           (solutionSP.UmlValue((ai_new-1)*param.Sdimmult+1)<JUcValue))
                            for k=1:4
                                sel='00';
                                binary = dec2bin(k-1);
                                sel(end-numel(binary)+1:end)=binary;
                                sel=reshape(str2num(char(num2cell(sel))),1,[]);
                                sel=sel+1;                  
                                i_self= (ai_new-1)*gridP.af.n*gridP.wf.n+(wi_m(sel(2))-1)*gridP.af.n+ai_m(sel(1));
                                AJ_J(i,i) = AJ_J(i,i)-param.m*solutionS_s.hl(ai_new)^1*(ai_m_w(sel(1))*wi_m_w(sel(2)))/singles;
                                AJ_JU(i,i_self)=AJ_JU(i,i_self)+param.m*solutionS_s.hl(ai_new)^1*(ai_m_w(sel(1))*wi_m_w(sel(2)))/singles;
                            end  

                        end
                        if always_marry||((solutionSP.JfcValue(ai, wi)<JUcValue)&&...
                           (solutionSP.UmcValue((ai_new-1)*param.Sdimmult+1)<JUcValue))
                            for k=1:4
                                sel='00';
                                binary = dec2bin(k-1);
                                sel(end-numel(binary)+1:end)=binary;
                                sel=reshape(str2num(char(num2cell(sel))),1,[]);
                                sel=sel+1;                  
                                i_self= (ai_new-1)*gridP.af.n*gridP.wf.n+(wi_m(sel(2))-1)*gridP.af.n+ai_m(sel(1));
                                AJ_J(i,i) = AJ_J(i,i)-param.m*solutionS_s.hc(ai_new)^1*(ai_m_w(sel(1))*wi_m_w(sel(2)))/singles;
                                AJ_JU(i,i_self)=AJ_JU(i,i_self)+param.m*solutionS_s.hc(ai_new)^1*(ai_m_w(sel(1))*wi_m_w(sel(2)))/singles;
                            end      
                        end

                        JOcValue=0;
                        for k=1:4
                            sel='00';
                            binary = dec2bin(k-1);
                            sel(end-numel(binary)+1:end)=binary;
                            sel=reshape(str2num(char(num2cell(sel))),1,[]);
                            sel=sel+1;                         
                            JOcValue=JOcValue+...
                                solutionSP.JOcValue(ai_m(sel(1)), wi_m(sel(2)), ai_new)*...
                                (ai_m_w(sel(1))*wi_m_w(sel(2)));                            
                        end
                        if always_marry||((solutionSP.JfcValue(ai, wi)<JOcValue)&&...
                           (solutionSP.OmlValue((ai_new-1)*param.Sdimmult+1)<JOcValue))
                            for k=1:4
                                sel='00';
                                binary = dec2bin(k-1);
                                sel(end-numel(binary)+1:end)=binary;
                                sel=reshape(str2num(char(num2cell(sel))),1,[]);
                                sel=sel+1;                  
                                i_self= (ai_new-1)*gridP.af.n*gridP.wf.n+(wi_m(sel(2))-1)*gridP.af.n+ai_m(sel(1));
                                AJ_J(i,i) = AJ_J(i,i)-param.m*solutionS_s.ol(ai_new)^1*(ai_m_w(sel(1))*wi_m_w(sel(2)))/singles;
                                AJ_JO(i,i_self)=AJ_JO(i,i_self)+param.m*solutionS_s.ol(ai_new)^1*(ai_m_w(sel(1))*wi_m_w(sel(2)))/singles;
                            end  

                        end
                        if always_marry||((solutionSP.JfcValue(ai, wi)<JOcValue)&&...
                           (solutionSP.OmcValue((ai_new-1)*param.Sdimmult+1)<JOcValue))
                            for k=1:4
                                sel='00';
                                binary = dec2bin(k-1);
                                sel(end-numel(binary)+1:end)=binary;
                                sel=reshape(str2num(char(num2cell(sel))),1,[]);
                                sel=sel+1;                  
                                i_self= (ai_new-1)*gridP.af.n*gridP.wf.n+(wi_m(sel(2))-1)*gridP.af.n+ai_m(sel(1));
                                AJ_J(i,i) = AJ_J(i,i)-param.m*solutionS_s.oc(ai_new)^1*(ai_m_w(sel(1))*wi_m_w(sel(2)))/singles;
                                AJ_JO(i,i_self)=AJ_JO(i,i_self)+param.m*solutionS_s.oc(ai_new)^1*(ai_m_w(sel(1))*wi_m_w(sel(2)))/singles;
                            end      
                        end
                    end                
                end
            end
        end
    end   
    AJ_J = sparse(AJ_J)+pg.CJ;
    AJ_U = pg.CJ_U;  
    AJ_O = sparse(AJ_O);
    AJ_JJ = sparse(AJ_JJ);
    AJ_JU = sparse(AJ_JU);
    AJ_JO = sparse(AJ_JO);
end
function [AJ_J, AJ_U, AJ_O,  AJ_JJ, AJ_UJ, AJ_OJ] = generateAJc_m(solutionSP,solutionS_s, pg, gridS)
    AJ_J = (zeros(pg.NS,pg.NS));
    AJ_O = zeros(size(pg.CJ_U));
    AJ_JJ = (zeros(pg.NS,param.NP));
    AJ_UJ = (zeros(pg.NS,param.NP/gridP.wf.n));
    AJ_OJ = (zeros(pg.NS,param.NP/gridP.wf.n));
    singles=1;%sum(solutionSP.h,'all')+sum(solutionSP.e,'all');

    for ai = 1:gridS.a.n
        %find weights for linear interpolation
        ai_m=mVai_m(ai,:);
        ai_m_w=mVai_m_w(ai,:);

        for wi = 1:gridS.w.n
            %find weights for linear interpolation
            wi_m=mVwi_m(wi,:);
            wi_m_w=mVwi_m_w(wi,:);

            i = (wi-1)*gridS.a.n + ai;
            if solutionSP.JmcValue(ai, wi) < solutionSP.OmcValue(ai)
                AJ_J(i,i)   = AJ_J(i,i)   - pg.lambda_o;
                AJ_O(i,ai)  = AJ_O(i,ai)  + pg.lambda_o;
            elseif tremble~=0
                t = tremble/exp( solutionSP.JmcValue(ai, wi) - solutionSP.OmcValue(ai) );
                AJ_J(i,i)   = AJ_J(i,i)   - t;
                AJ_O(i,ai)  = AJ_O(i,ai)  + t;
            end
            % for wi_new = 1:gridS.w.n
            %     i_new = (wi_new-1)*gridS.a.n + ai;
            %     if solutionSP.JmcValue(ai, wi)<solutionSP.JmcValue(ai, wi_new)
            %         AJ_J(i,i) = AJ_J(i,i)-pg.lambda*param.f_long(wi_new);
            %         AJ_J(i,i_new) = AJ_J(i,i_new)+pg.lambda*param.f_long(wi_new);
            %     end
            % end
            for wi_new = 1:gridP.wf.n
                %====
                %check marrige
                for ai_new = 1:gridP.af.n                                      
                    JJcValue=0;                    
                    for k=1:4
                        sel='00';
                        binary = dec2bin(k-1);
                        sel(end-numel(binary)+1:end)=binary;
                        sel=reshape(str2num(char(num2cell(sel))),1,[]);
                        sel=sel+1;                       
                        JJcValue=JJcValue+...
                            solutionSP.JJcValue(ai_new, wi_new, ai_m(sel(1)), wi_m(sel(2)))*...
                            (ai_m_w(sel(1))*wi_m_w(sel(2)));                        
                    end


                    if always_marry||((solutionSP.JmcValue(ai, wi)<JJcValue)&&...
                       (solutionSP.JflValue((ai_new-1)*param.Sdimmult+1, (wi_new-1)*param.Sdimmult+1)<JJcValue))
                        for k=1:4
                            sel='00';
                            binary = dec2bin(k-1);
                            sel(end-numel(binary)+1:end)=binary;
                            sel=reshape(str2num(char(num2cell(sel))),1,[]);
                            sel=sel+1;
                            i_self= (wi_m(sel(2))-1)*gridP.am.n*gridP.af.n*gridP.wf.n+...
                                (ai_m(sel(1))-1)*gridP.af.n*gridP.wf.n+...
                                (wi_new-1)*gridP.af.n+...
                                ai_new;
                            AJ_J(i,i) = AJ_J(i,i)-param.m*solutionS_s.el(ai_new,wi_new)^1*(ai_m_w(sel(1))*wi_m_w(sel(2)))/singles;
                            AJ_JJ(i,i_self)=AJ_JJ(i,i_self)+param.m*solutionS_s.el(ai_new,wi_new)^1*(ai_m_w(sel(1))*wi_m_w(sel(2)))/singles;
                        end                        
                    end
                    if always_marry||((solutionSP.JmcValue(ai, wi)<JJcValue)&&...
                       (solutionSP.JfcValue((ai_new-1)*param.Sdimmult+1, (wi_new-1)*param.Sdimmult+1)<JJcValue))
                        for k=1:4
                            sel='00';
                            binary = dec2bin(k-1);
                            sel(end-numel(binary)+1:end)=binary;
                            sel=reshape(str2num(char(num2cell(sel))),1,[]);
                            sel=sel+1;
                            i_self= (wi_m(sel(2))-1)*gridP.am.n*gridP.af.n*gridP.wf.n+...
                                (ai_m(sel(1))-1)*gridP.af.n*gridP.wf.n+...
                                (wi_new-1)*gridP.af.n+...
                                ai_new;
                            AJ_J(i,i) = AJ_J(i,i)-param.m*solutionS_s.ec(ai_new,wi_new)^1*(ai_m_w(sel(1))*wi_m_w(sel(2)))/singles;
                            AJ_JJ(i,i_self)=AJ_JJ(i,i_self)+param.m*solutionS_s.ec(ai_new,wi_new)^1*(ai_m_w(sel(1))*wi_m_w(sel(2)))/singles;
                        end   
                    end
                    if wi_new==1

                        UJcValue=0;
                        for k=1:4
                            sel='00';
                            binary = dec2bin(k-1);
                            sel(end-numel(binary)+1:end)=binary;
                            sel=reshape(str2num(char(num2cell(sel))),1,[]);
                            sel=sel+1;                          
                            UJcValue=UJcValue+...
                                solutionSP.UJcValue(ai_new, ai_m(sel(1)), wi_m(sel(2)))*...
                                (ai_m_w(sel(1))*wi_m_w(sel(2)));                            
                        end


                        if always_marry||((solutionSP.JmcValue(ai, wi)<UJcValue)&&...
                           (solutionSP.UflValue((ai_new-1)*param.Sdimmult+1)<UJcValue))
                            for k=1:4
                                sel='00';
                                binary = dec2bin(k-1);
                                sel(end-numel(binary)+1:end)=binary;
                                sel=reshape(str2num(char(num2cell(sel))),1,[]);
                                sel=sel+1;                  
                                i_self= (wi_m(sel(2))-1)*gridP.am.n*gridP.af.n+...
                                    (ai_m(sel(1))-1)*gridP.af.n+...
                                    ai_new;
                                AJ_J(i,i) = AJ_J(i,i)-param.m*solutionS_s.hl(ai_new)^1*(ai_m_w(sel(1))*wi_m_w(sel(2)))/singles;
                                AJ_UJ(i,i_self)=AJ_UJ(i,i_self)+param.m*solutionS_s.hl(ai_new)^1*(ai_m_w(sel(1))*wi_m_w(sel(2)))/singles;
                            end  

                        end
                        if always_marry||((solutionSP.JmcValue(ai, wi)<UJcValue)&&...
                           (solutionSP.UfcValue((ai_new-1)*param.Sdimmult+1)<UJcValue))
                            for k=1:4
                                sel='00';
                                binary = dec2bin(k-1);
                                sel(end-numel(binary)+1:end)=binary;
                                sel=reshape(str2num(char(num2cell(sel))),1,[]);
                                sel=sel+1;                  
                                i_self= (wi_m(sel(2))-1)*gridP.am.n*gridP.af.n+...
                                    (ai_m(sel(1))-1)*gridP.af.n+...
                                    ai_new;
                                AJ_J(i,i) = AJ_J(i,i)-param.m*solutionS_s.hc(ai_new)^1*(ai_m_w(sel(1))*wi_m_w(sel(2)))/singles;
                                AJ_UJ(i,i_self)=AJ_UJ(i,i_self)+param.m*solutionS_s.hc(ai_new)^1*(ai_m_w(sel(1))*wi_m_w(sel(2)))/singles;
                            end      
                        end

                        OJcValue=0;
                        for k=1:4
                            sel='00';
                            binary = dec2bin(k-1);
                            sel(end-numel(binary)+1:end)=binary;
                            sel=reshape(str2num(char(num2cell(sel))),1,[]);
                            sel=sel+1;                          
                            OJcValue=OJcValue+...
                                solutionSP.OJcValue(ai_new, ai_m(sel(1)), wi_m(sel(2)))*...
                                (ai_m_w(sel(1))*wi_m_w(sel(2)));                            
                        end
                        if always_marry||((solutionSP.JmcValue(ai, wi)<OJcValue)&&...
                           (solutionSP.OflValue((ai_new-1)*param.Sdimmult+1)<OJcValue))
                            for k=1:4
                                sel='00';
                                binary = dec2bin(k-1);
                                sel(end-numel(binary)+1:end)=binary;
                                sel=reshape(str2num(char(num2cell(sel))),1,[]);
                                sel=sel+1;                  
                                i_self= (wi_m(sel(2))-1)*gridP.am.n*gridP.af.n+...
                                    (ai_m(sel(1))-1)*gridP.af.n+...
                                    ai_new;
                                AJ_J(i,i) = AJ_J(i,i)-param.m*solutionS_s.ol(ai_new)^1*(ai_m_w(sel(1))*wi_m_w(sel(2)))/singles;
                                AJ_OJ(i,i_self)=AJ_OJ(i,i_self)+param.m*solutionS_s.ol(ai_new)^1*(ai_m_w(sel(1))*wi_m_w(sel(2)))/singles;
                            end  

                        end
                        if always_marry||((solutionSP.JmcValue(ai, wi)<OJcValue)&&...
                           (solutionSP.OfcValue((ai_new-1)*param.Sdimmult+1)<OJcValue))
                            for k=1:4
                                sel='00';
                                binary = dec2bin(k-1);
                                sel(end-numel(binary)+1:end)=binary;
                                sel=reshape(str2num(char(num2cell(sel))),1,[]);
                                sel=sel+1;                  
                                i_self= (wi_m(sel(2))-1)*gridP.am.n*gridP.af.n+...
                                    (ai_m(sel(1))-1)*gridP.af.n+...
                                    ai_new;
                                AJ_J(i,i) = AJ_J(i,i)-param.m*solutionS_s.oc(ai_new)^1*(ai_m_w(sel(1))*wi_m_w(sel(2)))/singles;
                                AJ_OJ(i,i_self)=AJ_OJ(i,i_self)+param.m*solutionS_s.oc(ai_new)^1*(ai_m_w(sel(1))*wi_m_w(sel(2)))/singles;
                            end      
                        end
                    end                
                end
            end
        end
    end   
    AJ_J = sparse(AJ_J)+pg.CJ;
    AJ_U = pg.CJ_U;  
    AJ_O = sparse(AJ_O);
    AJ_JJ = sparse(AJ_JJ);
    AJ_UJ = sparse(AJ_UJ);
    AJ_OJ = sparse(AJ_OJ);
end
function [AU_J, AU_U, AU_O, AU_UJ, AU_UU, AU_UO] = generateAUc_f(solutionSP,solutionS_s, pg, gridS)
    AU_J = (zeros(gridS.a.n,pg.NS));
    AU_U = (zeros(gridS.a.n,gridS.a.n));
    AU_O = zeros(size(AU_U));
    AU_UU = (zeros(gridS.a.n,gridP.af.n*gridP.am.n));
    AU_UO = (zeros(gridS.a.n,gridP.af.n*gridP.am.n));
    AU_UJ = (zeros(gridS.a.n,param.NP/gridP.wf.n));
    singles=1;%sum(solutionSP.h,'all')+sum(solutionSP.e,'all');
    for ai = 1:gridS.a.n
        %find weights for linear interpolation
        ai_m=fVai_m(ai,:);
        ai_m_w=fVai_m_w(ai,:);
        if solutionSP.UfcValue(ai) < solutionSP.OfcValue(ai)
            AU_U(ai,ai) = AU_U(ai,ai) - pg.lambda_o;
            AU_O(ai,ai) = AU_O(ai,ai) + pg.lambda_o;
        elseif tremble~=0
            t = tremble/exp( solutionSP.UfcValue(ai) - solutionSP.OfcValue(ai) );
            AU_U(ai,ai) = AU_U(ai,ai) - t;
            AU_O(ai,ai) = AU_O(ai,ai) + t;
        end
        for wi_new = 1:gridS.w.n                
            i_new = (wi_new-1)*gridS.a.n + ai;
            if true&&(solutionSP.UfcValue(ai)<solutionSP.JfcValue(ai, wi_new))%solutionSP.UfcValue(ai)<=(solutionSP.JfcValue(ai, wi_new)-param.crit)
                AU_U(ai,ai) = AU_U(ai,ai)-pg.lambda_u*param.f_long(wi_new);
                AU_J(ai,i_new) = AU_J(ai,i_new)+pg.lambda_u*param.f_long(wi_new);
            elseif tremble~=0
                t = tremble/exp(solutionSP.UfcValue(ai)-solutionSP.JfcValue(ai, wi_new));
                AU_U(ai,ai) = AU_U(ai,ai)-t*param.f_long(wi_new);
                AU_J(ai,i_new) = AU_J(ai,i_new)+t*param.f_long(wi_new);
            end
        end
        for wi_new = 1:gridP.wm.n     
            %check marrige
            for ai_new = 1:gridP.am.n                

                UJcValue=0;                
                for k=1:2
                    sel='0';
                    binary = dec2bin(k-1);
                    sel(end-numel(binary)+1:end)=binary;
                    sel=reshape(str2num(char(num2cell(sel))),1,[]);
                    sel=sel+1;                   
                    UJcValue=UJcValue+...
                        solutionSP.UJcValue(ai_m(sel(1)), ai_new, wi_new)*...
                        (ai_m_w(sel(1)));                    
                end

                if always_marry||((solutionSP.UfcValue(ai)<UJcValue)&&...
                   (solutionSP.JmlValue((ai_new-1)*param.Sdimmult+1, (wi_new-1)*param.Sdimmult+1)<UJcValue))    
                    for k=1:2
                        sel='0';
                        binary = dec2bin(k-1);
                        sel(end-numel(binary)+1:end)=binary;
                        sel=reshape(str2num(char(num2cell(sel))),1,[]);
                        sel=sel+1;
                        i_self = (wi_new-1)*gridP.af.n*gridP.am.n+(ai_new-1)*gridP.af.n+ai_m(sel(1));      
                        AU_U(ai,ai) = AU_U(ai,ai)-param.m*solutionS_s.el(ai_new,wi_new)^1*(ai_m_w(sel(1)))/singles;
                        AU_UJ(ai,i_self)=AU_UJ(ai,i_self)+param.m*solutionS_s.el(ai_new,wi_new)^1*(ai_m_w(sel(1)))/singles;
                    end                        
                end
                if always_marry||((solutionSP.UfcValue(ai)<UJcValue)&&...
                   (solutionSP.JmcValue((ai_new-1)*param.Sdimmult+1, (wi_new-1)*param.Sdimmult+1)<UJcValue))
                    for k=1:2
                        sel='0';
                        binary = dec2bin(k-1);
                        sel(end-numel(binary)+1:end)=binary;
                        sel=reshape(str2num(char(num2cell(sel))),1,[]);
                        sel=sel+1;
                        i_self = (wi_new-1)*gridP.af.n*gridP.am.n+(ai_new-1)*gridP.af.n+ai_m(sel(1));      
                        AU_U(ai,ai) = AU_U(ai,ai)-param.m*solutionS_s.ec(ai_new,wi_new)^1*(ai_m_w(sel(1)))/singles;
                        AU_UJ(ai,i_self)=AU_UJ(ai,i_self)+param.m*solutionS_s.ec(ai_new,wi_new)^1*(ai_m_w(sel(1)))/singles;
                    end   
                end
                if wi_new==1

                    UUcValue=0;
                    for k=1:2
                        sel='0';
                        binary = dec2bin(k-1);
                        sel(end-numel(binary)+1:end)=binary;
                        sel=reshape(str2num(char(num2cell(sel))),1,[]);
                        sel=sel+1;                      
                        UUcValue=UUcValue+...
                            solutionSP.UUcValue(ai_m(sel(1)), ai_new)*...
                            (ai_m_w(sel(1)));                        
                    end


                    if always_marry||((solutionSP.UfcValue(ai)<UUcValue)&&...
                       (solutionSP.UmlValue((ai_new-1)*param.Sdimmult+1)<UUcValue))
                        for k=1:2
                            sel='0';
                            binary = dec2bin(k-1);
                            sel(end-numel(binary)+1:end)=binary;
                            sel=reshape(str2num(char(num2cell(sel))),1,[]);
                            sel=sel+1;                  
                            i_self= (ai_new-1)*gridP.af.n+ai_m(sel(1));
                            AU_U(ai,ai) = AU_U(ai,ai)-param.m*solutionS_s.hl(ai_new)^1*(ai_m_w(sel(1)))/singles;
                            AU_UU(ai,i_self)=AU_UU(ai,i_self)+param.m*solutionS_s.hl(ai_new)^1*(ai_m_w(sel(1)))/singles;
                        end                          
                    end
                    if always_marry||((solutionSP.UfcValue(ai)<UUcValue)&&...
                       (solutionSP.UmcValue((ai_new-1)*param.Sdimmult+1)<UUcValue))
                        for k=1:2
                            sel='0';
                            binary = dec2bin(k-1);
                            sel(end-numel(binary)+1:end)=binary;
                            sel=reshape(str2num(char(num2cell(sel))),1,[]);
                            sel=sel+1;                  
                            i_self= (ai_new-1)*gridP.af.n+ai_m(sel(1));
                            AU_U(ai,ai) = AU_U(ai,ai)-param.m*solutionS_s.hc(ai_new)^1*(ai_m_w(sel(1)))/singles;
                            AU_UU(ai,i_self)=AU_UU(ai,i_self)+param.m*solutionS_s.hc(ai_new)^1*(ai_m_w(sel(1)))/singles;
                        end     
                    end

                    UOcValue=0;
                    for k=1:2
                        sel='0';
                        binary = dec2bin(k-1);
                        sel(end-numel(binary)+1:end)=binary;
                        sel=reshape(str2num(char(num2cell(sel))),1,[]);
                        sel=sel+1;                      
                        UOcValue=UOcValue+...
                            solutionSP.UOcValue(ai_m(sel(1)), ai_new)*...
                            (ai_m_w(sel(1)));                        
                    end
                    if always_marry||((solutionSP.UfcValue(ai)<UOcValue)&&...
                       (solutionSP.OmlValue((ai_new-1)*param.Sdimmult+1)<UOcValue))
                        for k=1:2
                            sel='0';
                            binary = dec2bin(k-1);
                            sel(end-numel(binary)+1:end)=binary;
                            sel=reshape(str2num(char(num2cell(sel))),1,[]);
                            sel=sel+1;                  
                            i_self= (ai_new-1)*gridP.af.n+ai_m(sel(1));
                            AU_U(ai,ai) = AU_U(ai,ai)-param.m*solutionS_s.ol(ai_new)^1*(ai_m_w(sel(1)))/singles;
                            AU_UO(ai,i_self)=AU_UO(ai,i_self)+param.m*solutionS_s.ol(ai_new)^1*(ai_m_w(sel(1)))/singles;
                        end                          
                    end

                    if always_marry||((solutionSP.UfcValue(ai)<UOcValue)&&...
                       (solutionSP.OmcValue((ai_new-1)*param.Sdimmult+1)<UOcValue))
                        for k=1:2
                            sel='0';
                            binary = dec2bin(k-1);
                            sel(end-numel(binary)+1:end)=binary;
                            sel=reshape(str2num(char(num2cell(sel))),1,[]);
                            sel=sel+1;                  
                            i_self= (ai_new-1)*gridP.af.n+ai_m(sel(1));
                            AU_U(ai,ai) = AU_U(ai,ai)-param.m*solutionS_s.oc(ai_new)^1*(ai_m_w(sel(1)))/singles;
                            AU_UO(ai,i_self)=AU_UO(ai,i_self)+param.m*solutionS_s.oc(ai_new)^1*(ai_m_w(sel(1)))/singles;
                        end     
                    end
                end

            end
        end
    end 
    AU_J = sparse(AU_J);
    AU_U = sparse(AU_U) + pg.CU;
    AU_O = sparse(AU_O);
    AU_UU = sparse(AU_UU);
    AU_UO = sparse(AU_UO);
    AU_UJ = sparse(AU_UJ);
end
function [AU_J, AU_U, AU_O, AU_JU, AU_UU, AU_OU] = generateAUc_m(solutionSP,solutionS_s, pg, gridS)
    AU_J = (zeros(gridS.a.n,pg.NS));
    AU_U = (zeros(gridS.a.n,gridS.a.n));
    AU_O = zeros(size(AU_U));
    AU_UU = (zeros(gridS.a.n,gridP.af.n*gridP.am.n));
    AU_OU = (zeros(gridS.a.n,gridP.af.n*gridP.am.n));
    AU_JU = (zeros(gridS.a.n,param.NP/gridP.wm.n));
    singles=1;%sum(solutionSP.h,'all')+sum(solutionSP.e,'all');
    for ai = 1:gridS.a.n
        %find weights for linear interpolation
        ai_m=mVai_m(ai,:);
        ai_m_w=mVai_m_w(ai,:);
        if solutionSP.UmcValue(ai) < solutionSP.OmcValue(ai)
            AU_U(ai,ai) = AU_U(ai,ai) - pg.lambda_o;
            AU_O(ai,ai) = AU_O(ai,ai) + pg.lambda_o;
        elseif tremble~=0
            t = tremble/exp( solutionSP.UmcValue(ai) - solutionSP.OmcValue(ai) );
            AU_U(ai,ai) = AU_U(ai,ai) - t;
            AU_O(ai,ai) = AU_O(ai,ai) + t;
        end
        for wi_new = 1:gridS.w.n                
            i_new = (wi_new-1)*gridS.a.n + ai;
            if true&&(solutionSP.UmcValue(ai)<solutionSP.JmcValue(ai, wi_new))%solutionSP.UmcValue(ai)<=(solutionSP.JmcValue(ai, wi_new)-param.crit)
                AU_U(ai,ai) = AU_U(ai,ai)-pg.lambda_u*param.f_long(wi_new);
                AU_J(ai,i_new) = AU_J(ai,i_new)+pg.lambda_u*param.f_long(wi_new);
            elseif tremble~=0
                t = tremble/exp(solutionSP.UmcValue(ai)-solutionSP.JmcValue(ai, wi_new));
                AU_U(ai,ai) = AU_U(ai,ai)-t*param.f_long(wi_new);
                AU_J(ai,i_new) = AU_J(ai,i_new)+t*param.f_long(wi_new);
            end
        end
        for wi_new = 1:gridP.wf.n     %was wm
            %check marrige
            for ai_new = 1:gridP.am.n                

                JUcValue=0;                
                for k=1:2
                    sel='0';
                    binary = dec2bin(k-1);
                    sel(end-numel(binary)+1:end)=binary;
                    sel=reshape(str2num(char(num2cell(sel))),1,[]);
                    sel=sel+1;               
                    JUcValue=JUcValue+...
                        solutionSP.JUcValue(ai_new, wi_new, ai_m(sel(1)))*...
                        (ai_m_w(sel(1)));                    
                end

                if always_marry||((solutionSP.UmcValue(ai)<JUcValue)&&...
                   (solutionSP.JflValue((ai_new-1)*param.Sdimmult+1, (wi_new-1)*param.Sdimmult+1)<JUcValue))
                    for k=1:2
                        sel='0';
                        binary = dec2bin(k-1);
                        sel(end-numel(binary)+1:end)=binary;
                        sel=reshape(str2num(char(num2cell(sel))),1,[]);
                        sel=sel+1;             
                        i_self= (ai_m(sel(1))-1)*gridP.af.n*gridP.wf.n+...
                                (wi_new-1)*gridP.af.n+...
                                ai_new;
                        AU_U(ai,ai) = AU_U(ai,ai)-param.m*solutionS_s.el(ai_new,wi_new)^1*(ai_m_w(sel(1)))/singles;
                        AU_JU(ai,i_self)=AU_JU(ai,i_self)+param.m*solutionS_s.el(ai_new,wi_new)^1*(ai_m_w(sel(1)))/singles;
                    end                        
                end
                if always_marry||((solutionSP.UmcValue(ai)<JUcValue)&&...
                   (solutionSP.JfcValue((ai_new-1)*param.Sdimmult+1, (wi_new-1)*param.Sdimmult+1)<JUcValue))
                    for k=1:2
                        sel='0';
                        binary = dec2bin(k-1);
                        sel(end-numel(binary)+1:end)=binary;
                        sel=reshape(str2num(char(num2cell(sel))),1,[]);
                        sel=sel+1;
                        i_self= (ai_m(sel(1))-1)*gridP.af.n*gridP.wf.n+...
                                (wi_new-1)*gridP.af.n+...
                                ai_new;      
                        AU_U(ai,ai) = AU_U(ai,ai)-param.m*solutionS_s.ec(ai_new,wi_new)^1*(ai_m_w(sel(1)))/singles;
                        AU_JU(ai,i_self)=AU_JU(ai,i_self)+param.m*solutionS_s.ec(ai_new,wi_new)^1*(ai_m_w(sel(1)))/singles;
                    end   
                end
                if wi_new==1

                    UUcValue=0;
                    for k=1:2
                        sel='0';
                        binary = dec2bin(k-1);
                        sel(end-numel(binary)+1:end)=binary;
                        sel=reshape(str2num(char(num2cell(sel))),1,[]);
                        sel=sel+1;                      
                        UUcValue=UUcValue+...
                            solutionSP.UUcValue(ai_m(sel(1)), ai_new)*...
                            (ai_m_w(sel(1)));                        
                    end


                    if always_marry||((solutionSP.UmcValue(ai)<UUcValue)&&...
                       (solutionSP.UflValue((ai_new-1)*param.Sdimmult+1)<UUcValue))    
                        for k=1:2
                            sel='0';
                            binary = dec2bin(k-1);
                            sel(end-numel(binary)+1:end)=binary;
                            sel=reshape(str2num(char(num2cell(sel))),1,[]);
                            sel=sel+1;                  
                            i_self= (ai_m(sel(1))-1)*gridP.af.n+ai_new;
                            AU_U(ai,ai) = AU_U(ai,ai)-param.m*solutionS_s.hl(ai_new)^1*(ai_m_w(sel(1)))/singles;
                            AU_UU(ai,i_self)=AU_UU(ai,i_self)+param.m*solutionS_s.hl(ai_new)^1*(ai_m_w(sel(1)))/singles;
                        end                          
                    end
                    if always_marry||((solutionSP.UmcValue(ai)<UUcValue)&&...
                       (solutionSP.UfcValue((ai_new-1)*param.Sdimmult+1)<UUcValue))  
                        for k=1:2
                            sel='0';
                            binary = dec2bin(k-1);
                            sel(end-numel(binary)+1:end)=binary;
                            sel=reshape(str2num(char(num2cell(sel))),1,[]);
                            sel=sel+1;                  
                            i_self= (ai_m(sel(1))-1)*gridP.af.n+ai_new;
                            AU_U(ai,ai) = AU_U(ai,ai)-param.m*solutionS_s.hc(ai_new)^1*(ai_m_w(sel(1)))/singles;
                            AU_UU(ai,i_self)=AU_UU(ai,i_self)+param.m*solutionS_s.hc(ai_new)^1*(ai_m_w(sel(1)))/singles;
                        end     
                    end

                    OUcValue=0;
                    for k=1:2
                        sel='0';
                        binary = dec2bin(k-1);
                        sel(end-numel(binary)+1:end)=binary;
                        sel=reshape(str2num(char(num2cell(sel))),1,[]);
                        sel=sel+1;                      
                        OUcValue=OUcValue+...
                            solutionSP.OUcValue(ai_m(sel(1)), ai_new)*...
                            (ai_m_w(sel(1)));                        
                    end
                    if always_marry||((solutionSP.UmcValue(ai)<OUcValue)&&...
                       (solutionSP.OflValue((ai_new-1)*param.Sdimmult+1)<OUcValue))    
                        for k=1:2
                            sel='0';
                            binary = dec2bin(k-1);
                            sel(end-numel(binary)+1:end)=binary;
                            sel=reshape(str2num(char(num2cell(sel))),1,[]);
                            sel=sel+1;                  
                            i_self= (ai_m(sel(1))-1)*gridP.af.n+ai_new;
                            AU_U(ai,ai) = AU_U(ai,ai)-param.m*solutionS_s.ol(ai_new)^1*(ai_m_w(sel(1)))/singles;
                            AU_OU(ai,i_self)=AU_OU(ai,i_self)+param.m*solutionS_s.ol(ai_new)^1*(ai_m_w(sel(1)))/singles;
                        end                          
                    end
                    if always_marry||((solutionSP.UmcValue(ai)<OUcValue)&&...
                       (solutionSP.OfcValue((ai_new-1)*param.Sdimmult+1)<OUcValue))  
                        for k=1:2
                            sel='0';
                            binary = dec2bin(k-1);
                            sel(end-numel(binary)+1:end)=binary;
                            sel=reshape(str2num(char(num2cell(sel))),1,[]);
                            sel=sel+1;                  
                            i_self= (ai_m(sel(1))-1)*gridP.af.n+ai_new;
                            AU_U(ai,ai) = AU_U(ai,ai)-param.m*solutionS_s.oc(ai_new)^1*(ai_m_w(sel(1)))/singles;
                            AU_OU(ai,i_self)=AU_OU(ai,i_self)+param.m*solutionS_s.oc(ai_new)^1*(ai_m_w(sel(1)))/singles;
                        end     
                    end
                end

            end
        end
    end 
    AU_J = sparse(AU_J);
    AU_U = sparse(AU_U) + pg.CU;
    AU_O = sparse(AU_O);
    AU_UU = sparse(AU_UU);
    AU_OU = sparse(AU_OU);
    AU_JU = sparse(AU_JU);
end
function [AO_O, AO_J, AO_OJ, AO_OU] = generateAOc_f(solutionSP,solutionS_s, pg, gridS)    
    AO_O  = (zeros(gridS.a.n,gridS.a.n));
    AO_J  = (zeros(gridS.a.n,pg.NS));              % NEW: O -> J (single employed with kids)
    AO_OU = (zeros(gridS.a.n,gridP.af.n*gridP.am.n));
    AO_OJ = (zeros(gridS.a.n,param.NP/gridP.wf.n));
    singles=1; % scaling as in your code

    for ai = 1:gridS.a.n
        ai_m   = fVai_m(ai,:);
        ai_m_w = fVai_m_w(ai,:);

        % ---- O -> J over the single wage grid (replaces O -> U) ----
        for wi_new = 1:gridS.w.n
            i_new = (wi_new-1)*gridS.a.n + ai;   % destination: single J(state) with kids
            if solutionSP.OfcValue(ai) < solutionSP.JfcValue(ai, wi_new)
                AO_O(ai,ai)  = AO_O(ai,ai)  - param.pf.lambda_oE * param.f_long(wi_new);
                AO_J(ai,i_new)= AO_J(ai,i_new)+ param.pf.lambda_oE * param.f_long(wi_new);
            elseif tremble~=0
                t = tremble/exp( solutionSP.OfcValue(ai) - solutionSP.JfcValue(ai, wi_new) );
                AO_O(ai,ai)  = AO_O(ai,ai)  - t * param.f_long(wi_new);
                AO_J(ai,i_new)= AO_J(ai,i_new)+ t * param.f_long(wi_new);
            end
        end

        for wi_new = 1:gridP.wm.n     
            %check marrige
            for ai_new = 1:gridP.am.n                

                OJcValue=0;                
                for k=1:2
                    sel='0';
                    binary = dec2bin(k-1);
                    sel(end-numel(binary)+1:end)=binary;
                    sel=reshape(str2num(char(num2cell(sel))),1,[]);
                    sel=sel+1;                   
                    OJcValue=OJcValue+...
                        solutionSP.OJcValue(ai_m(sel(1)), ai_new, wi_new)*...
                        (ai_m_w(sel(1)));                    
                end

                if always_marry||((solutionSP.OfcValue(ai)<OJcValue)&&...
                   (solutionSP.JmlValue((ai_new-1)*param.Sdimmult+1, (wi_new-1)*param.Sdimmult+1)<OJcValue))   
                    for k=1:2
                        sel='0';
                        binary = dec2bin(k-1);
                        sel(end-numel(binary)+1:end)=binary;
                        sel=reshape(str2num(char(num2cell(sel))),1,[]);
                        sel=sel+1;
                        i_self = (wi_new-1)*gridP.af.n*gridP.am.n+(ai_new-1)*gridP.af.n+ai_m(sel(1));      
                        AO_O(ai,ai) = AO_O(ai,ai)-param.m*solutionS_s.el(ai_new,wi_new)^1*(ai_m_w(sel(1)))/singles;
                        AO_OJ(ai,i_self)=AO_OJ(ai,i_self)+param.m*solutionS_s.el(ai_new,wi_new)^1*(ai_m_w(sel(1)))/singles;
                    end                        
                end
                if always_marry||((solutionSP.OfcValue(ai)<OJcValue)&&...
                   (solutionSP.JmcValue((ai_new-1)*param.Sdimmult+1, (wi_new-1)*param.Sdimmult+1)<OJcValue))    
                    for k=1:2
                        sel='0';
                        binary = dec2bin(k-1);
                        sel(end-numel(binary)+1:end)=binary;
                        sel=reshape(str2num(char(num2cell(sel))),1,[]);
                        sel=sel+1;
                        i_self = (wi_new-1)*gridP.af.n*gridP.am.n+(ai_new-1)*gridP.af.n+ai_m(sel(1));      
                        AO_O(ai,ai) = AO_O(ai,ai)-param.m*solutionS_s.ec(ai_new,wi_new)^1*(ai_m_w(sel(1)))/singles;
                        AO_OJ(ai,i_self)=AO_OJ(ai,i_self)+param.m*solutionS_s.ec(ai_new,wi_new)^1*(ai_m_w(sel(1)))/singles;
                    end   
                end
                if wi_new==1

                    OUcValue=0;
                    for k=1:2
                        sel='0';
                        binary = dec2bin(k-1);
                        sel(end-numel(binary)+1:end)=binary;
                        sel=reshape(str2num(char(num2cell(sel))),1,[]);
                        sel=sel+1;

                        OUcValue=OUcValue+...
                            solutionSP.OUcValue(ai_m(sel(1)), ai_new)*...
                            (ai_m_w(sel(1)));                        
                    end


                    if always_marry||((solutionSP.OfcValue(ai)<OUcValue)&&...
                       (solutionSP.UmlValue((ai_new-1)*param.Sdimmult+1)<OUcValue))    
                        for k=1:2
                            sel='0';
                            binary = dec2bin(k-1);
                            sel(end-numel(binary)+1:end)=binary;
                            sel=reshape(str2num(char(num2cell(sel))),1,[]);
                            sel=sel+1;                  
                            i_self= (ai_new-1)*gridP.af.n+ai_m(sel(1));
                            AO_O(ai,ai) = AO_O(ai,ai)-param.m*solutionS_s.hl(ai_new)^1*(ai_m_w(sel(1)))/singles;
                            AO_OU(ai,i_self)=AO_OU(ai,i_self)+param.m*solutionS_s.hl(ai_new)^1*(ai_m_w(sel(1)))/singles;
                        end                          
                    end
                    if always_marry||((solutionSP.OfcValue(ai)<OUcValue)&&...
                       (solutionSP.UmcValue((ai_new-1)*param.Sdimmult+1)<OUcValue))    
                        for k=1:2
                            sel='0';
                            binary = dec2bin(k-1);
                            sel(end-numel(binary)+1:end)=binary;
                            sel=reshape(str2num(char(num2cell(sel))),1,[]);
                            sel=sel+1;                  
                            i_self= (ai_new-1)*gridP.af.n+ai_m(sel(1));
                            AO_O(ai,ai) = AO_O(ai,ai)-param.m*solutionS_s.hc(ai_new)^1*(ai_m_w(sel(1)))/singles;
                            AO_OU(ai,i_self)=AO_OU(ai,i_self)+param.m*solutionS_s.hc(ai_new)^1*(ai_m_w(sel(1)))/singles;
                        end     
                    end
                end
            end
        end
    end     
    AO_O = sparse(AO_O) + pg.CO;
    AO_J = sparse(AO_J);
    AO_OU = sparse(AO_OU);
    AO_OJ = sparse(AO_OJ);
end

function [AO_O, AO_J, AO_JO, AO_UO] = generateAOc_m(solutionSP,solutionS_s, pg, gridS)    
    AO_O  = (zeros(gridS.a.n,gridS.a.n));
    AO_J  = (zeros(gridS.a.n,pg.NS));              % NEW: O -> J (single employed with kids)
    AO_UO = (zeros(gridS.a.n,gridP.af.n*gridP.am.n));
    AO_JO = (zeros(gridS.a.n,param.NP/gridP.wm.n));
    singles=1;

    for ai = 1:gridS.a.n
        ai_m   = mVai_m(ai,:);
        ai_m_w = mVai_m_w(ai,:);

        % ---- O -> J over the single wage grid (replaces O -> U) ----
        for wi_new = 1:gridS.w.n
            i_new = (wi_new-1)*gridS.a.n + ai;   % destination: single J(state) with kids
            if solutionSP.OmcValue(ai) < solutionSP.JmcValue(ai, wi_new)
                AO_O(ai,ai)  = AO_O(ai,ai)  - param.pm.lambda_oE * param.f_long(wi_new);
                AO_J(ai,i_new)= AO_J(ai,i_new)+ param.pm.lambda_oE * param.f_long(wi_new);
            elseif tremble~=0
                t = tremble/exp( solutionSP.OmcValue(ai) - solutionSP.JmcValue(ai, wi_new) );
                AO_O(ai,ai)  = AO_O(ai,ai)  - t * param.f_long(wi_new);
                AO_J(ai,i_new)= AO_J(ai,i_new)+ t * param.f_long(wi_new);
            end
        end
        for wi_new = 1:gridP.wf.n % was wm     
            %check marrige
            for ai_new = 1:gridP.am.n                

                JOcValue=0;                
                for k=1:2
                    sel='0';
                    binary = dec2bin(k-1);
                    sel(end-numel(binary)+1:end)=binary;
                    sel=reshape(str2num(char(num2cell(sel))),1,[]);
                    sel=sel+1;

                    JOcValue=JOcValue+...
                        solutionSP.JOcValue(ai_new, wi_new, ai_m(sel(1)))*...
                        (ai_m_w(sel(1)));                    
                end

                if always_marry||((solutionSP.OmcValue(ai)<JOcValue)&&...
                   (solutionSP.JflValue((ai_new-1)*param.Sdimmult+1, (wi_new-1)*param.Sdimmult+1)<JOcValue))    
                    for k=1:2
                        sel='0';
                        binary = dec2bin(k-1);
                        sel(end-numel(binary)+1:end)=binary;
                        sel=reshape(str2num(char(num2cell(sel))),1,[]);
                        sel=sel+1;             
                        i_self= (ai_m(sel(1))-1)*gridP.af.n*gridP.wf.n+...
                                (wi_new-1)*gridP.af.n+...
                                ai_new;
                        AO_O(ai,ai) = AO_O(ai,ai)-param.m*solutionS_s.el(ai_new,wi_new)^1*(ai_m_w(sel(1)))/singles;
                        AO_JO(ai,i_self)=AO_JO(ai,i_self)+param.m*solutionS_s.el(ai_new,wi_new)^1*(ai_m_w(sel(1)))/singles;
                    end                        
                end
                if always_marry||((solutionSP.OmcValue(ai)<JOcValue)&&...
                   (solutionSP.JfcValue((ai_new-1)*param.Sdimmult+1, (wi_new-1)*param.Sdimmult+1)<JOcValue))    
                    for k=1:2
                        sel='0';
                        binary = dec2bin(k-1);
                        sel(end-numel(binary)+1:end)=binary;
                        sel=reshape(str2num(char(num2cell(sel))),1,[]);
                        sel=sel+1;
                        i_self= (ai_m(sel(1))-1)*gridP.af.n*gridP.wf.n+...
                                (wi_new-1)*gridP.af.n+...
                                ai_new;      
                        AO_O(ai,ai) = AO_O(ai,ai)-param.m*solutionS_s.ec(ai_new,wi_new)^1*(ai_m_w(sel(1)))/singles;
                        AO_JO(ai,i_self)=AO_JO(ai,i_self)+param.m*solutionS_s.ec(ai_new,wi_new)^1*(ai_m_w(sel(1)))/singles;
                    end   
                end
                if wi_new==1

                    UOcValue=0;
                    for k=1:2
                        sel='0';
                        binary = dec2bin(k-1);
                        sel(end-numel(binary)+1:end)=binary;
                        sel=reshape(str2num(char(num2cell(sel))),1,[]);
                        sel=sel+1;

                        UOcValue=UOcValue+...
                            solutionSP.UOcValue(ai_m(sel(1)), ai_new)*...
                            (ai_m_w(sel(1)));                        
                    end


                    if always_marry||((solutionSP.OmcValue(ai)<UOcValue)&&...
                       (solutionSP.UflValue((ai_new-1)*param.Sdimmult+1)<UOcValue))   
                        for k=1:2
                            sel='0';
                            binary = dec2bin(k-1);
                            sel(end-numel(binary)+1:end)=binary;
                            sel=reshape(str2num(char(num2cell(sel))),1,[]);
                            sel=sel+1;                  
                            i_self= (ai_m(sel(1))-1)*gridP.af.n+ai_new;
                            AO_O(ai,ai) = AO_O(ai,ai)-param.m*solutionS_s.hl(ai_new)^1*(ai_m_w(sel(1)))/singles;
                            AO_UO(ai,i_self)=AO_UO(ai,i_self)+param.m*solutionS_s.hl(ai_new)^1*(ai_m_w(sel(1)))/singles;
                        end                          
                    end
                    if always_marry||((solutionSP.OmcValue(ai)<UOcValue)&&...
                       (solutionSP.UfcValue((ai_new-1)*param.Sdimmult+1)<UOcValue))    
                        for k=1:2
                            sel='0';
                            binary = dec2bin(k-1);
                            sel(end-numel(binary)+1:end)=binary;
                            sel=reshape(str2num(char(num2cell(sel))),1,[]);
                            sel=sel+1;                  
                            i_self= (ai_m(sel(1))-1)*gridP.af.n+ai_new;
                            AO_O(ai,ai) = AO_O(ai,ai)-param.m*solutionS_s.hc(ai_new)^1*(ai_m_w(sel(1)))/singles;
                            AO_UO(ai,i_self)=AO_UO(ai,i_self)+param.m*solutionS_s.hc(ai_new)^1*(ai_m_w(sel(1)))/singles;
                        end     
                    end
                end

            end
        end
    end     
    AO_O = sparse(AO_O) + pg.CO;
    AO_J = sparse(AO_J);
    AO_UO = sparse(AO_UO);
    AO_JO = sparse(AO_JO);
end

function [AJJ_JJ, AJJ_JU, AJJ_UJ, AJJ_JO, AJJ_OJ] = generateAJJl(solutionSP, param, gridP)
    % Dimensions
    na  = gridP.af.n;
    nwf = gridP.wf.n;
    nam = gridP.am.n;
    nwm = gridP.wm.n;
    NP  = param.NP;

    % Prealloc with explicit sizes
    AJJ_JJ = zeros(NP, NP);                         % JJ -> JJ
    AJJ_JU = param.CJJ_JU;                          % JJ -> JU (constant block, provided)
    AJJ_UJ = param.CJJ_UJ;                          % JJ -> UJ (constant block, provided)
    AJJ_JO = zeros(NP, NP / nwm);                   % JJ -> JO (drop wm)
    AJJ_OJ = zeros(NP, NP / nwf);                   % JJ -> OJ (drop wf)

    for afi = 1:na
        for wfi = 1:nwf
            for ami = 1:nam
                for wmi = 1:nwm
                    i = (wmi-1)*nam*na*nwf + (ami-1)*na*nwf + (wfi-1)*na + afi;

                    % --- JJ -> OJ (female to O): drops wf dimension
                    i_OJ = (wmi-1)*nam*na + (ami-1)*na + afi;  % (wm, am, af)
                    if solutionSP.JJlValue(afi,wfi,ami,wmi) < solutionSP.OJlValue(afi,ami,wmi)
                        AJJ_JJ(i,i)   = AJJ_JJ(i,i)   - param.pf.lambda_o;
                        AJJ_OJ(i,i_OJ)= AJJ_OJ(i,i_OJ)+ param.pf.lambda_o;
                    elseif tremble ~= 0
                        t = tremble / exp(solutionSP.JJlValue(afi,wfi,ami,wmi) - solutionSP.OJlValue(afi,ami,wmi));
                        AJJ_JJ(i,i)   = AJJ_JJ(i,i)   - t;
                        AJJ_OJ(i,i_OJ)= AJJ_OJ(i,i_OJ)+ t;
                    end

                    % --- JJ -> JO (male to O): drops wm dimension
                    i_JO = (ami-1)*na*nwf + (wfi-1)*na + afi;  % (am, wf, af)
                    if solutionSP.JJlValue(afi,wfi,ami,wmi) < solutionSP.JOlValue(afi,wfi,ami)
                        AJJ_JJ(i,i)   = AJJ_JJ(i,i)   - param.pm.lambda_o;
                        AJJ_JO(i,i_JO)= AJJ_JO(i,i_JO)+ param.pm.lambda_o;
                    elseif tremble ~= 0
                        t = tremble / exp(solutionSP.JJlValue(afi,wfi,ami,wmi) - solutionSP.JOlValue(afi,wfi,ami));
                        AJJ_JJ(i,i)   = AJJ_JJ(i,i)   - t;
                        AJJ_JO(i,i_JO)= AJJ_JO(i,i_JO)+ t;
                    end
                end
            end
        end
    end

    % Finalize sparsity / add constant diagonal block
    AJJ_JJ = sparse(AJJ_JJ) + param.CJJ;
    AJJ_JO = sparse(AJJ_JO);
    AJJ_OJ = sparse(AJJ_OJ);
end


function [AJU_JU, AJU_JJ, AJU_UU, AJU_OU, AJU_JO] = generateAJUl(solutionSP, param, gridP)
    AJU_JU = zeros(param.NP/gridP.wm.n, param.NP/gridP.wm.n);
    AJU_JJ = zeros(param.NP/gridP.wm.n, param.NP);
    AJU_OU = zeros(size(param.CJU_UU));
    AJU_JO = zeros(size(AJU_JU));

    for afi = 1:gridP.af.n
        for wfi = 1:gridP.wf.n
            for ami = 1:gridP.am.n
                i = (ami-1)*gridP.af.n*gridP.wf.n + (wfi-1)*gridP.af.n + afi;

                % --- JU -> OU (female J->O): uses pf.lambda_o
                i_new = (ami-1)*gridP.af.n + afi;    % (af, am) ordering
                if solutionSP.JUlValue(afi,wfi,ami) < solutionSP.OUlValue(afi,ami)
                    AJU_JU(i,i)     = AJU_JU(i,i)     - param.pf.lambda_o;
                    AJU_OU(i,i_new) = AJU_OU(i,i_new) + param.pf.lambda_o;
                elseif tremble~=0
                    t = tremble/exp(solutionSP.JUlValue(afi,wfi,ami) - solutionSP.OUlValue(afi,ami));
                    AJU_JU(i,i)     = AJU_JU(i,i)     - t;
                    AJU_OU(i,i_new) = AJU_OU(i,i_new) + t;
                end

                % --- JU -> JO (male U->O): uses pm.lambda_o
                if solutionSP.JUlValue(afi,wfi,ami) < solutionSP.JOlValue(afi,wfi,ami)
                    AJU_JU(i,i) = AJU_JU(i,i) - param.pm.lambda_o;
                    AJU_JO(i,i) = AJU_JO(i,i) + param.pm.lambda_o;
                elseif tremble~=0
                    t = tremble/exp(solutionSP.JUlValue(afi,wfi,ami) - solutionSP.JOlValue(afi,wfi,ami));
                    AJU_JU(i,i) = AJU_JU(i,i) - t;
                    AJU_JO(i,i) = AJU_JO(i,i) + t;
                end

                % --- JU -> JJ (male job search over wm)
                for wmi_new = 1:gridP.wm.n
                    i_new = (wmi_new-1)*gridP.am.n*gridP.af.n*gridP.wf.n + ...
                            (ami-1)*gridP.af.n*gridP.wf.n + (wfi-1)*gridP.af.n + afi;
                    if solutionSP.JUlValue(afi,wfi,ami) < solutionSP.JJlValue(afi,wfi,ami,wmi_new)
                        AJU_JU(i,i)     = AJU_JU(i,i)     - param.pm.lambda_u * param.f(wmi_new);
                        AJU_JJ(i,i_new) = AJU_JJ(i,i_new) + param.pm.lambda_u * param.f(wmi_new);
                    elseif tremble~=0
                        t = tremble/exp(solutionSP.JUlValue(afi,wfi,ami) - solutionSP.JJlValue(afi,wfi,ami,wmi_new));
                        AJU_JU(i,i)     = AJU_JU(i,i)     - t * param.f(wmi_new);
                        AJU_JJ(i,i_new) = AJU_JJ(i,i_new) + t * param.f(wmi_new);
                    end
                end
            end
        end
    end

    AJU_JU = sparse(AJU_JU) + param.CJU;
    AJU_JJ = sparse(AJU_JJ);
    AJU_UU = param.CJU_UU;
    AJU_JO = sparse(AJU_JO);
    AJU_OU = sparse(AJU_OU);
end

    function [AJO_JO, AJO_UO, AJO_OO, AJO_JJ] = generateAJOl(solutionSP, param, gridP)
    AJO_JO = (zeros(param.NP/gridP.wm.n,param.NP/gridP.wm.n));
    AJO_OO = zeros(size(param.CJU_UU));
    AJO_JJ = (zeros(param.NP/gridP.wm.n, param.NP));

    for afi = 1:gridP.af.n
        for wfi = 1:gridP.wf.n
            for ami = 1:gridP.am.n
                i = (ami-1)*gridP.af.n*gridP.wf.n+(wfi-1)*gridP.af.n+afi;

                % JO -> OO (female quits to OLF): unchanged
                i_new_OO = (ami-1)*gridP.af.n+afi;
                if solutionSP.JOlValue(afi,wfi,ami) < solutionSP.OOlValue(afi,ami)
                    AJO_JO(i,i)   = AJO_JO(i,i)   - param.pf.lambda_o;
                    AJO_OO(i,i_new_OO) = AJO_OO(i,i_new_OO) + param.pf.lambda_o;
                elseif tremble~=0
                    t = tremble/exp( solutionSP.JOlValue(afi,wfi,ami) - solutionSP.OOlValue(afi,ami) );
                    AJO_JO(i,i)   = AJO_JO(i,i)   - t;
                    AJO_OO(i,i_new_OO) = AJO_OO(i,i_new_OO) + t;
                end

                % REPLACEMENT: JO -> JJ via male OLF → E (over wm)
                for wmi_new = 1:gridP.wm.n
                    i_new_JJ = (wmi_new-1)*gridP.am.n*gridP.af.n*gridP.wf.n + ...
                               (ami-1)*gridP.af.n*gridP.wf.n + (wfi-1)*gridP.af.n + afi;
                    if solutionSP.JOlValue(afi,wfi,ami) < solutionSP.JJlValue(afi,wfi,ami,wmi_new)
                        AJO_JO(i,i)     = AJO_JO(i,i)     - param.pm.lambda_oE * param.f(wmi_new);
                        AJO_JJ(i,i_new_JJ) = AJO_JJ(i,i_new_JJ) + param.pm.lambda_oE * param.f(wmi_new);
                    elseif tremble~=0
                        t = tremble/exp( solutionSP.JOlValue(afi,wfi,ami) - solutionSP.JJlValue(afi,wfi,ami,wmi_new) );
                        AJO_JO(i,i)     = AJO_JO(i,i)     - t * param.f(wmi_new);
                        AJO_JJ(i,i_new_JJ) = AJO_JJ(i,i_new_JJ) + t * param.f(wmi_new);
                    end
                end
            end
        end
    end
    AJO_JO = sparse(AJO_JO)+ param.CJU;
    AJO_UO = param.CJU_UU;
    AJO_OO = sparse(AJO_OO);
    AJO_JJ = sparse(AJO_JJ);
end

function [AUJ_UJ, AUJ_JJ, AUJ_UU, AUJ_OJ, AUJ_UO] = generateAUJl(solutionSP, param, gridP)
    AUJ_UJ = zeros(param.NP/gridP.wf.n, param.NP/gridP.wf.n);
    AUJ_JJ = zeros(param.NP/gridP.wf.n, param.NP);
    AUJ_OJ = zeros(size(AUJ_UJ));
    AUJ_UO = zeros(size(param.CUJ_UU));

    for afi = 1:gridP.af.n
        for ami = 1:gridP.am.n
            for wmi = 1:gridP.wm.n
                i = (wmi-1)*gridP.af.n*gridP.am.n + (ami-1)*gridP.af.n + afi;

                % --- UJ -> UO (male U->O): pm.lambda_o  
                i_new = (ami-1)*gridP.af.n + afi;   % (af, am)
                if solutionSP.UJlValue(afi,ami,wmi) < solutionSP.UOlValue(afi,ami)
                    AUJ_UJ(i,i)     = AUJ_UJ(i,i)     - param.pm.lambda_o;   
                    AUJ_UO(i,i_new) = AUJ_UO(i,i_new) + param.pm.lambda_o;
                elseif tremble~=0
                    t = tremble/exp(solutionSP.UJlValue(afi,ami,wmi) - solutionSP.UOlValue(afi,ami));
                    AUJ_UJ(i,i)     = AUJ_UJ(i,i)     - t;
                    AUJ_UO(i,i_new) = AUJ_UO(i,i_new) + t;
                end

                % --- UJ -> OJ (female J->O): pf.lambda_o 
                if solutionSP.UJlValue(afi,ami,wmi) < solutionSP.OJlValue(afi,ami,wmi)
                    AUJ_UJ(i,i) = AUJ_UJ(i,i) - param.pf.lambda_o;           
                    AUJ_OJ(i,i) = AUJ_OJ(i,i) + param.pf.lambda_o;
                elseif tremble~=0
                    t = tremble/exp(solutionSP.UJlValue(afi,ami,wmi) - solutionSP.OJlValue(afi,ami,wmi));
                    AUJ_UJ(i,i) = AUJ_UJ(i,i) - t;
                    AUJ_OJ(i,i) = AUJ_OJ(i,i) + t;
                end

                % --- UJ -> JJ (female job search over wf)
                for wfi_new = 1:gridP.wf.n
                    i_new = (wmi-1)*gridP.am.n*gridP.af.n*gridP.wf.n + ...
                            (ami-1)*gridP.af.n*gridP.wf.n + (wfi_new-1)*gridP.af.n + afi;
                    if solutionSP.UJlValue(afi,ami,wmi) < solutionSP.JJlValue(afi,wfi_new,ami,wmi)
                        AUJ_UJ(i,i)     = AUJ_UJ(i,i)     - param.pf.lambda_u * param.f(wfi_new);
                        AUJ_JJ(i,i_new) = AUJ_JJ(i,i_new) + param.pf.lambda_u * param.f(wfi_new);
                    elseif tremble~=0
                        t = tremble/exp(solutionSP.UJlValue(afi,ami,wmi) - solutionSP.JJlValue(afi,wfi_new,ami,wmi));
                        AUJ_UJ(i,i)     = AUJ_UJ(i,i)     - t * param.f(wfi_new);
                        AUJ_JJ(i,i_new) = AUJ_JJ(i,i_new) + t * param.f(wfi_new);
                    end
                end
            end
        end
    end

    AUJ_UJ = sparse(AUJ_UJ) + param.CUJ;
    AUJ_JJ = sparse(AUJ_JJ);
    AUJ_UU = param.CUJ_UU;
    AUJ_OJ = sparse(AUJ_OJ);
    AUJ_UO = sparse(AUJ_UO);
end

function [AOJ_OJ, AOJ_OU, AOJ_OO, AOJ_JJ] = generateAOJl(solutionSP, param, gridP)
    AOJ_OJ = (zeros(param.NP/gridP.wf.n,param.NP/gridP.wf.n));    
    AOJ_OO = zeros(size(param.CUJ_UU));
    % DEST now JJ (rows OJ = NP/wf, cols JJ = NP)
    AOJ_JJ = (zeros(param.NP/gridP.wf.n, param.NP));

    for afi = 1:gridP.af.n
        for ami = 1:gridP.am.n
            for wmi = 1:gridP.wm.n
                i = (wmi-1)*gridP.af.n*gridP.am.n+(ami-1)*gridP.af.n+afi;                     

                % OJ -> OO (male quits to OLF): unchanged
                i_new_OO = (ami-1)*gridP.af.n+afi;        
                if solutionSP.OJlValue(afi,ami,wmi) < solutionSP.OOlValue(afi,ami)
                    AOJ_OJ(i,i)   = AOJ_OJ(i,i)   - param.pm.lambda_o;
                    AOJ_OO(i,i_new_OO) = AOJ_OO(i,i_new_OO) + param.pm.lambda_o;
                elseif tremble~=0
                    t = tremble/exp( solutionSP.OJlValue(afi,ami,wmi) - solutionSP.OOlValue(afi,ami) );
                    AOJ_OJ(i,i)   = AOJ_OJ(i,i)   - t;
                    AOJ_OO(i,i_new_OO) = AOJ_OO(i,i_new_OO) + t;
                end
                
                % REPLACEMENT: OJ -> JJ via female OLF → E (over wf)
                for wfi_new = 1:gridP.wf.n
                    i_new_JJ = (wmi-1)*gridP.am.n*gridP.af.n*gridP.wf.n + ...
                               (ami-1)*gridP.af.n*gridP.wf.n + (wfi_new-1)*gridP.af.n + afi;
                    if solutionSP.OJlValue(afi,ami,wmi) < solutionSP.JJlValue(afi,wfi_new,ami,wmi)
                        AOJ_OJ(i,i)       = AOJ_OJ(i,i)       - param.pf.lambda_oE * param.f(wfi_new);
                        AOJ_JJ(i,i_new_JJ) = AOJ_JJ(i,i_new_JJ) + param.pf.lambda_oE * param.f(wfi_new);
                    elseif tremble~=0
                        t = tremble/exp( solutionSP.OJlValue(afi,ami,wmi) - solutionSP.JJlValue(afi,wfi_new,ami,wmi) );
                        AOJ_OJ(i,i)       = AOJ_OJ(i,i)       - t * param.f(wfi_new);
                        AOJ_JJ(i,i_new_JJ) = AOJ_JJ(i,i_new_JJ) + t * param.f(wfi_new);
                    end
                end
            end
        end        
    end 
    AOJ_OJ = sparse(AOJ_OJ)+ param.CUJ;
    AOJ_OU = param.CUJ_UU;
    AOJ_OO = sparse(AOJ_OO);
    AOJ_JJ = sparse(AOJ_JJ);
end

function [AUU_UU, AUU_UJ, AUU_JU, AUU_UO, AUU_OU] = generateAUUl(solutionSP, param, gridP)
    AUU_UU = zeros(gridP.af.n*gridP.am.n, gridP.af.n*gridP.am.n);
    AUU_UJ = zeros(gridP.af.n*gridP.am.n, param.NP/gridP.wf.n);
    AUU_JU = zeros(gridP.af.n*gridP.am.n, param.NP/gridP.wm.n);
    AUU_UO = zeros(size(AUU_UU));
    AUU_OU = zeros(size(AUU_UU));
    for afi = 1:gridP.af.n
        for ami = 1:gridP.am.n
            i = (ami-1)*gridP.af.n + afi;

            % UU -> UO (male U→O): female hazard  **FIXED**
            if solutionSP.UUlValue(afi,ami) < solutionSP.UOlValue(afi,ami)
                AUU_UU(i,i) = AUU_UU(i,i) - param.pm.lambda_o;   % was pm.lambda_o
                AUU_UO(i,i) = AUU_UO(i,i) + param.pm.lambda_o;
            elseif tremble~=0
                t = tremble/exp(solutionSP.UUlValue(afi,ami) - solutionSP.UOlValue(afi,ami));
                AUU_UU(i,i) = AUU_UU(i,i) - t;
                AUU_UO(i,i) = AUU_UO(i,i) + t;
            end

            % UU -> OU (female U→O): male hazard  **FIXED**
            if solutionSP.UUlValue(afi,ami) < solutionSP.OUlValue(afi,ami)
                AUU_UU(i,i) = AUU_UU(i,i) - param.pf.lambda_o;   % was pf.lambda_o
                AUU_OU(i,i) = AUU_OU(i,i) + param.pf.lambda_o;
            elseif tremble~=0
                t = tremble/exp(solutionSP.UUlValue(afi,ami) - solutionSP.OUlValue(afi,ami));
                AUU_UU(i,i) = AUU_UU(i,i) - t;
                AUU_OU(i,i) = AUU_OU(i,i) + t;
            end

            % UU -> JU (female job search over wf)
            for wfi_new = 1:gridP.wf.n
                i_new = (ami-1)*gridP.af.n*gridP.wf.n + (wfi_new-1)*gridP.af.n + afi;
                if solutionSP.UUlValue(afi,ami) < solutionSP.JUlValue(afi,wfi_new,ami)
                    AUU_UU(i,i)     = AUU_UU(i,i)     - param.pf.lambda_u * param.f(wfi_new);
                    AUU_JU(i,i_new) = AUU_JU(i,i_new) + param.pf.lambda_u * param.f(wfi_new);
                elseif tremble~=0
                    t = tremble/exp(solutionSP.UUlValue(afi,ami) - solutionSP.JUlValue(afi,wfi_new,ami));
                    AUU_UU(i,i)     = AUU_UU(i,i)     - t * param.f(wfi_new);
                    AUU_JU(i,i_new) = AUU_JU(i,i_new) + t * param.f(wfi_new);
                end
            end

            % UU -> UJ (male job search over wm)
            for wmi_new = 1:gridP.wm.n
                i_new = (wmi_new-1)*gridP.af.n*gridP.am.n + (ami-1)*gridP.af.n + afi;
                if solutionSP.UUlValue(afi,ami) < solutionSP.UJlValue(afi,ami,wmi_new)
                    AUU_UU(i,i)     = AUU_UU(i,i)     - param.pm.lambda_u * param.f(wmi_new);
                    AUU_UJ(i,i_new) = AUU_UJ(i,i_new) + param.pm.lambda_u * param.f(wmi_new);
                elseif tremble~=0
                    t = tremble/exp(solutionSP.UUlValue(afi,ami) - solutionSP.UJlValue(afi,ami,wmi_new));
                    AUU_UU(i,i)     = AUU_UU(i,i)     - t * param.f(wmi_new);
                    AUU_UJ(i,i_new) = AUU_UJ(i,i_new) + t * param.f(wmi_new);
                end
            end
        end
    end
    AUU_UU = sparse(AUU_UU) + param.CUU;
    AUU_UJ = sparse(AUU_UJ);
    AUU_JU = sparse(AUU_JU);
    AUU_UO = sparse(AUU_UO);
    AUU_OU = sparse(AUU_OU);
end

function [AOU_OU, AOU_OJ, AOU_OO, AOU_JU] = generateAOUl(solutionSP, param, gridP)
    AOU_OU = (zeros(gridP.af.n*gridP.am.n,gridP.af.n*gridP.am.n));
    AOU_OJ = (zeros(gridP.af.n*gridP.am.n,param.NP/gridP.wf.n));    
    AOU_OO = zeros(size(AOU_OU));
    % DEST now JU (rows OU = af*am, cols JU = NP/wm)
    AOU_JU = (zeros(gridP.af.n*gridP.am.n, param.NP/gridP.wm.n));

    for afi = 1:gridP.af.n
        for ami = 1:gridP.am.n
            i = (ami-1)*gridP.af.n+afi; 

            % OU -> OO (male U→O): unchanged
            if solutionSP.OUlValue(afi,ami) < solutionSP.OOlValue(afi,ami)
                AOU_OU(i,i) = AOU_OU(i,i) - param.pm.lambda_o;
                AOU_OO(i,i) = AOU_OO(i,i) + param.pm.lambda_o;
            elseif tremble~=0
                t = tremble/exp( solutionSP.OUlValue(afi,ami) - solutionSP.OOlValue(afi,ami) );
                AOU_OU(i,i) = AOU_OU(i,i) - t;
                AOU_OO(i,i) = AOU_OO(i,i) + t;
            end

            % REPLACEMENT: OU -> JU via female OLF → E (over wf)
            for wfi_new = 1:gridP.wf.n
                i_new_JU = (ami-1)*gridP.af.n*gridP.wf.n + (wfi_new-1)*gridP.af.n + afi;
                if solutionSP.OUlValue(afi,ami) < solutionSP.JUlValue(afi,wfi_new,ami)
                    AOU_OU(i,i)     = AOU_OU(i,i)     - param.pf.lambda_oE * param.f(wfi_new);
                    AOU_JU(i,i_new_JU) = AOU_JU(i,i_new_JU) + param.pf.lambda_oE * param.f(wfi_new);
                elseif tremble~=0
                    t = tremble/exp( solutionSP.OUlValue(afi,ami) - solutionSP.JUlValue(afi,wfi_new,ami) );
                    AOU_OU(i,i)     = AOU_OU(i,i)     - t * param.f(wfi_new);
                    AOU_JU(i,i_new_JU) = AOU_JU(i,i_new_JU) + t * param.f(wfi_new);
                end
            end

            % male U→J (OU → OJ): unchanged
            for wmi_new = 1:gridP.wm.n
                i_new_OJ = (wmi_new-1)*gridP.af.n*gridP.am.n + (ami-1)*gridP.af.n + afi;
                if true&&(solutionSP.OUlValue(afi,ami)<solutionSP.OJlValue(afi,ami,wmi_new))
                    AOU_OU(i,i) = AOU_OU(i,i)-param.pm.lambda_u*param.f(wmi_new);
                    AOU_OJ(i,i_new_OJ) = AOU_OJ(i,i_new_OJ)+param.pm.lambda_u*param.f(wmi_new);
                elseif tremble~=0
                    t = tremble/exp(solutionSP.OUlValue(afi,ami)-solutionSP.OJlValue(afi,ami,wmi_new));
                    AOU_OU(i,i) = AOU_OU(i,i)-t*param.f(wmi_new);
                    AOU_OJ(i,i_new_OJ) = AOU_OJ(i,i_new_OJ)+t*param.f(wmi_new);
                end
            end
        end        
    end 
    AOU_OU = sparse(AOU_OU) + param.CUU;
    AOU_OJ = sparse(AOU_OJ);
    AOU_OO = sparse(AOU_OO);
    AOU_JU = sparse(AOU_JU);
end

function [AUO_UO, AUO_JO, AUO_OO, AUO_UJ] = generateAUOl(solutionSP, param, gridP)
    AUO_UO = (zeros(gridP.af.n*gridP.am.n,gridP.af.n*gridP.am.n));    
    AUO_JO = (zeros(gridP.af.n*gridP.am.n,param.NP/gridP.wf.n));
    AUO_OO = zeros(size(AUO_UO));
    % DEST now UJ (rows UO = af*am, cols UJ = NP/wf)
    AUO_UJ = (zeros(gridP.af.n*gridP.am.n, param.NP/gridP.wf.n));

    for afi = 1:gridP.af.n
        for ami = 1:gridP.am.n
            i = (ami-1)*gridP.af.n+afi; 

            % UO -> OO (female U→O): unchanged
            if solutionSP.UOlValue(afi,ami) < solutionSP.OOlValue(afi,ami)
                AUO_UO(i,i) = AUO_UO(i,i) - param.pf.lambda_o;
                AUO_OO(i,i) = AUO_OO(i,i) + param.pf.lambda_o;
            elseif tremble~=0
                t = tremble/exp( solutionSP.UOlValue(afi,ami) - solutionSP.OOlValue(afi,ami) );
                AUO_UO(i,i) = AUO_UO(i,i) - t;
                AUO_OO(i,i) = AUO_OO(i,i) + t;
            end
            
            % REPLACEMENT: UO -> UJ via male OLF → E (over wm)
            for wmi_new = 1:gridP.wm.n
                i_new_UJ = (wmi_new-1)*gridP.af.n*gridP.am.n + (ami-1)*gridP.af.n + afi;
                if solutionSP.UOlValue(afi,ami) < solutionSP.UJlValue(afi,ami,wmi_new)
                    AUO_UO(i,i)      = AUO_UO(i,i)      - param.pm.lambda_oE * param.f(wmi_new);
                    AUO_UJ(i,i_new_UJ) = AUO_UJ(i,i_new_UJ) + param.pm.lambda_oE * param.f(wmi_new);
                elseif tremble~=0
                    t = tremble/exp( solutionSP.UOlValue(afi,ami) - solutionSP.UJlValue(afi,ami,wmi_new) );
                    AUO_UO(i,i)      = AUO_UO(i,i)      - t * param.f(wmi_new);
                    AUO_UJ(i,i_new_UJ) = AUO_UJ(i,i_new_UJ) + t * param.f(wmi_new);
                end
            end

            % female U→J (UO → JO): unchanged
            for wfi_new = 1:gridP.wf.n
                i_new_JO = (ami-1)*gridP.af.n*gridP.wf.n+(wfi_new-1)*gridP.af.n+afi;
                if true&&(solutionSP.UOlValue(afi,ami)<solutionSP.JOlValue(afi,wfi_new,ami))
                    AUO_UO(i,i) = AUO_UO(i,i)-param.pf.lambda_u*param.f(wfi_new);
                    AUO_JO(i,i_new_JO) = AUO_JO(i,i_new_JO)+param.pf.lambda_u*param.f(wfi_new);
                elseif tremble~=0
                    t = tremble/exp(solutionSP.UOlValue(afi,ami)-solutionSP.JOlValue(afi,wfi_new,ami));
                    AUO_UO(i,i) = AUO_UO(i,i)-t*param.f(wfi_new);
                    AUO_JO(i,i_new_JO) = AUO_JO(i,i_new_JO)+t*param.f(wfi_new);
                end
            end
        end        
    end 
    AUO_UO = sparse(AUO_UO) + param.CUU;
    AUO_JO = sparse(AUO_JO);
    AUO_OO = sparse(AUO_OO);
    AUO_UJ = sparse(AUO_UJ);
end

function [AOO_OO, AOO_OJ, AOO_JO] = generateAOOl(solutionSP, param, gridP)
    AOO_OO = (zeros(gridP.af.n*gridP.am.n,gridP.af.n*gridP.am.n));        
    % DESTS now OJ and JO (rows OO = af*am)
    AOO_OJ = (zeros(gridP.af.n*gridP.am.n, param.NP/gridP.wf.n)); % male to E -> OJ
    AOO_JO = (zeros(gridP.af.n*gridP.am.n, param.NP/gridP.wm.n)); % female to E -> JO

    for afi = 1:gridP.af.n
        for ami = 1:gridP.am.n
            i = (ami-1)*gridP.af.n+afi; 

            % REPLACEMENT: OO -> JO via female OLF → E (over wf)
            for wfi_new = 1:gridP.wf.n
                i_new_JO = (ami-1)*gridP.af.n*gridP.wf.n + (wfi_new-1)*gridP.af.n + afi;
                if solutionSP.OOlValue(afi,ami) < solutionSP.JOlValue(afi,wfi_new,ami)
                    AOO_OO(i,i)   = AOO_OO(i,i)   - param.pf.lambda_oE * param.f(wfi_new);
                    AOO_JO(i,i_new_JO) = AOO_JO(i,i_new_JO) + param.pf.lambda_oE * param.f(wfi_new);
                elseif tremble~=0
                    t = tremble/exp( solutionSP.OOlValue(afi,ami) - solutionSP.JOlValue(afi,wfi_new,ami) );
                    AOO_OO(i,i)   = AOO_OO(i,i)   - t * param.f(wfi_new);
                    AOO_JO(i,i_new_JO) = AOO_JO(i,i_new_JO) + t * param.f(wfi_new);
                end
            end

            % REPLACEMENT: OO -> OJ via male OLF → E (over wm)
            for wmi_new = 1:gridP.wm.n
                i_new_OJ = (wmi_new-1)*gridP.af.n*gridP.am.n + (ami-1)*gridP.af.n + afi;
                if solutionSP.OOlValue(afi,ami) < solutionSP.OJlValue(afi,ami,wmi_new)
                    AOO_OO(i,i)   = AOO_OO(i,i)   - param.pm.lambda_oE * param.f(wmi_new);
                    AOO_OJ(i,i_new_OJ) = AOO_OJ(i,i_new_OJ) + param.pm.lambda_oE * param.f(wmi_new);
                elseif tremble~=0
                    t = tremble/exp( solutionSP.OOlValue(afi,ami) - solutionSP.OJlValue(afi,ami,wmi_new) );
                    AOO_OO(i,i)   = AOO_OO(i,i)   - t * param.f(wmi_new);
                    AOO_OJ(i,i_new_OJ) = AOO_OJ(i,i_new_OJ) + t * param.f(wmi_new);
                end
            end
        end        
    end 
    AOO_OO = sparse(AOO_OO) + param.CUU;    
    AOO_OJ = sparse(AOO_OJ);
    AOO_JO = sparse(AOO_JO);
end

function [AJJ_JJ, AJJ_JU, AJJ_UJ, AJJ_JO, AJJ_OJ] = generateAJJc(solutionSP, param, gridP)
    % Dimensions
    na  = gridP.af.n;
    nwf = gridP.wf.n;
    nam = gridP.am.n;
    nwm = gridP.wm.n;
    NP  = param.NP;

    % Prealloc with explicit sizes
    AJJ_JJ = zeros(NP, NP);                         % JJc -> JJc
    AJJ_JU = param.CJJ_JU;                          % JJc -> JUc (constant)
    AJJ_UJ = param.CJJ_UJ;                          % JJc -> UJc (constant)
    AJJ_JO = zeros(NP, NP / nwm);                   % JJc -> JOc
    AJJ_OJ = zeros(NP, NP / nwf);                   % JJc -> OJc

    for afi = 1:na
        for wfi = 1:nwf
            for ami = 1:nam
                for wmi = 1:nwm
                    % Row index for JJc (wm, am, wf, af)
                    i = (wmi-1)*nam*na*nwf + (ami-1)*na*nwf + (wfi-1)*na + afi;

                    % --- JJc -> OJc (female to O)
                    i_OJ = (wmi-1)*nam*na + (ami-1)*na + afi;  % (wm, am, af)
                    if solutionSP.JJcValue(afi,wfi,ami,wmi) < solutionSP.OJcValue(afi,ami,wmi)
                        AJJ_JJ(i,i)   = AJJ_JJ(i,i)   - param.pf.lambda_o;
                        AJJ_OJ(i,i_OJ)= AJJ_OJ(i,i_OJ)+ param.pf.lambda_o;
                    elseif tremble ~= 0
                        t = tremble / exp(solutionSP.JJcValue(afi,wfi,ami,wmi) - solutionSP.OJcValue(afi,ami,wmi));
                        AJJ_JJ(i,i)   = AJJ_JJ(i,i)   - t;
                        AJJ_OJ(i,i_OJ)= AJJ_OJ(i,i_OJ)+ t;
                    end

                    % --- JJc -> JOc (male to O)
                    i_JO = (ami-1)*na*nwf + (wfi-1)*na + afi;  % (am, wf, af)
                    if solutionSP.JJcValue(afi,wfi,ami,wmi) < solutionSP.JOcValue(afi,wfi,ami)
                        AJJ_JJ(i,i)   = AJJ_JJ(i,i)   - param.pm.lambda_o;
                        AJJ_JO(i,i_JO)= AJJ_JO(i,i_JO)+ param.pm.lambda_o;
                    elseif tremble ~= 0
                        t = tremble / exp(solutionSP.JJcValue(afi,wfi,ami,wmi) - solutionSP.JOcValue(afi,wfi,ami));
                        AJJ_JJ(i,i)   = AJJ_JJ(i,i)   - t;
                        AJJ_JO(i,i_JO)= AJJ_JO(i,i_JO)+ t;
                    end
                end
            end
        end
    end

    % Finalize
    AJJ_JJ = sparse(AJJ_JJ) + param.CJJ;
    AJJ_JO = sparse(AJJ_JO);
    AJJ_OJ = sparse(AJJ_OJ);
end

function [AJU_JU, AJU_JJ, AJU_UU, AJU_OU, AJU_JO] = generateAJUc(solutionSP, param, gridP)
    % Dimensions
    na  = gridP.af.n;  nwf = gridP.wf.n;  nam = gridP.am.n;

    % Blocks
    AJU_JU = zeros(param.NP/gridP.wm.n, param.NP/gridP.wm.n);
    AJU_JJ = zeros(param.NP/gridP.wm.n, param.NP);
    AJU_OU = zeros(size(param.CJU_UU));      % (JU -> OU): (af,wf,am) x (af,am)
    AJU_JO = zeros(size(AJU_JU));            % (JU -> JO): (af,wf,am) x (af,wf,am)

    for afi = 1:na
        for wfi = 1:nwf
            for ami = 1:nam
                % Indices
                i     = (ami-1)*na*nwf + (wfi-1)*na + afi;    % (af,wf,am)
                i_OU  = (ami-1)*na + afi;                     % (af,am)

                % JU^c -> OU^c (female hazard)
                if solutionSP.JUcValue(afi,wfi,ami) < solutionSP.OUcValue(afi,ami)
                    AJU_JU(i,i)   = AJU_JU(i,i)   - param.pf.lambda_o;
                    AJU_OU(i,i_OU)= AJU_OU(i,i_OU)+ param.pf.lambda_o;
                elseif tremble ~= 0
                    t = tremble/exp(solutionSP.JUcValue(afi,wfi,ami) - solutionSP.OUcValue(afi,ami));
                    AJU_JU(i,i)   = AJU_JU(i,i)   - t;
                    AJU_OU(i,i_OU)= AJU_OU(i,i_OU)+ t;
                end

                % JU^c -> JO^c (male hazard)
                if solutionSP.JUcValue(afi,wfi,ami) < solutionSP.JOcValue(afi,wfi,ami)
                    AJU_JU(i,i) = AJU_JU(i,i) - param.pm.lambda_o;
                    AJU_JO(i,i) = AJU_JO(i,i) + param.pm.lambda_o;
                elseif tremble ~= 0
                    t = tremble/exp(solutionSP.JUcValue(afi,wfi,ami) - solutionSP.JOcValue(afi,wfi,ami));
                    AJU_JU(i,i) = AJU_JU(i,i) - t;
                    AJU_JO(i,i) = AJU_JO(i,i) + t;
                end

                % JU^c -> JJ^c (male job arrival over w_m)
                for wmi_new = 1:gridP.wm.n
                    jJJ = (wmi_new-1)*nam*na*nwf + (ami-1)*na*nwf + (wfi-1)*na + afi; % (af,wf,am,wm)
                    if solutionSP.JUcValue(afi,wfi,ami) < solutionSP.JJcValue(afi,wfi,ami,wmi_new)
                        AJU_JU(i,i)   = AJU_JU(i,i)   - param.pm.lambda_u * param.f(wmi_new);
                        AJU_JJ(i,jJJ) = AJU_JJ(i,jJJ) + param.pm.lambda_u * param.f(wmi_new);
                    elseif tremble ~= 0
                        t = tremble/exp(solutionSP.JUcValue(afi,wfi,ami) - solutionSP.JJcValue(afi,wfi,ami,wmi_new));
                        AJU_JU(i,i)   = AJU_JU(i,i)   - t * param.f(wmi_new);
                        AJU_JJ(i,jJJ) = AJU_JJ(i,jJJ) + t * param.f(wmi_new);
                    end
                end
            end
        end
    end

    AJU_JU = sparse(AJU_JU) + param.CJU;
    AJU_JJ = sparse(AJU_JJ);
    AJU_UU = param.CJU_UU;
    AJU_JO = sparse(AJU_JO);
    AJU_OU = sparse(AJU_OU);
end

function [AJO_JO, AJO_UO, AJO_OO, AJO_JJ] = generateAJOc(solutionSP, param, gridP)
    AJO_JO = (zeros(param.NP/gridP.wm.n,param.NP/gridP.wm.n));    
    AJO_OO = zeros(size(param.CJU_UU));
    AJO_JJ = (zeros(param.NP/gridP.wm.n,param.NP));

    for afi = 1:gridP.af.n
        for wfi = 1:gridP.wf.n
            for ami = 1:gridP.am.n
                i = (ami-1)*gridP.af.n*gridP.wf.n+(wfi-1)*gridP.af.n+afi;

                % JOc -> OOc (female quits to OLF): unchanged
                i_new = (ami-1)*gridP.af.n+afi;        
                if solutionSP.JOcValue(afi, wfi, ami) < solutionSP.OOcValue(afi, ami)
                    AJO_JO(i,i)     = AJO_JO(i,i)     - param.pf.lambda_o;
                    AJO_OO(i,i_new) = AJO_OO(i,i_new) + param.pf.lambda_o;
                elseif tremble~=0
                    t = tremble/exp( solutionSP.JOcValue(afi, wfi, ami) - solutionSP.OOcValue(afi, ami) );
                    AJO_JO(i,i)     = AJO_JO(i,i)     - t;
                    AJO_OO(i,i_new) = AJO_OO(i,i_new) + t;
                end

                % NEW: JOc -> JJc via male OLF -> E (over wmi)
                for wmi_new = 1:gridP.wm.n
                    i_new_JJ = (wmi_new-1)*gridP.am.n*gridP.af.n*gridP.wf.n + ...
                               (ami-1)*gridP.af.n*gridP.wf.n + (wfi-1)*gridP.af.n + afi;
                    if solutionSP.JOcValue(afi, wfi, ami) < solutionSP.JJcValue(afi, wfi, ami, wmi_new)
                        AJO_JO(i,i)         = AJO_JO(i,i)         - param.pm.lambda_oE * param.f(wmi_new);
                        AJO_JJ(i,i_new_JJ)  = AJO_JJ(i,i_new_JJ)  + param.pm.lambda_oE * param.f(wmi_new);
                    elseif tremble~=0
                        t = tremble/exp( solutionSP.JOcValue(afi, wfi, ami) - solutionSP.JJcValue(afi, wfi, ami, wmi_new) );
                        AJO_JO(i,i)         = AJO_JO(i,i)         - t * param.f(wmi_new);
                        AJO_JJ(i,i_new_JJ)  = AJO_JJ(i,i_new_JJ)  + t * param.f(wmi_new);
                    end
                end
            end
        end        
    end 
    AJO_JO = sparse(AJO_JO)+ param.CJU;    
    AJO_UO = param.CJU_UU;
    AJO_OO = sparse(AJO_OO);
    AJO_JJ = sparse(AJO_JJ);
end

function [AUJ_UJ, AUJ_JJ, AUJ_UU, AUJ_OJ, AUJ_UO] = generateAUJc(solutionSP, param, gridP)
    % Dimensions
    na  = gridP.af.n;  nam = gridP.am.n;  nwm = gridP.wm.n;

    % Blocks
    AUJ_UJ = zeros(param.NP/gridP.wf.n, param.NP/gridP.wf.n); % (UJ -> UJ)
    AUJ_JJ = zeros(param.NP/gridP.wf.n, param.NP);            % (UJ -> JJ)
    AUJ_OJ = zeros(size(AUJ_UJ));                              % (UJ -> OJ)
    AUJ_UO = zeros(size(param.CUJ_UU));                        % (UJ -> UO): (af,am,wm) x (af,am)

    for afi = 1:na
        for ami = 1:nam
            for wmi = 1:nwm
                i    = (wmi-1)*na*nam + (ami-1)*na + afi; % (af,am,wm)
                i_UO = (ami-1)*na + afi;                  % (af,am)

                % UJ^c -> UO^c (male hazard)
                if solutionSP.UJcValue(afi,ami,wmi) < solutionSP.UOcValue(afi,ami)
                    AUJ_UJ(i,i)   = AUJ_UJ(i,i)   - param.pm.lambda_o;
                    AUJ_UO(i,i_UO)= AUJ_UO(i,i_UO)+ param.pm.lambda_o;
                elseif tremble ~= 0
                    t = tremble/exp(solutionSP.UJcValue(afi,ami,wmi) - solutionSP.UOcValue(afi,ami));
                    AUJ_UJ(i,i)   = AUJ_UJ(i,i)   - t;
                    AUJ_UO(i,i_UO)= AUJ_UO(i,i_UO)+ t;
                end

                % UJ^c -> OJ^c (female hazard)
                if solutionSP.UJcValue(afi,ami,wmi) < solutionSP.OJcValue(afi,ami,wmi)
                    AUJ_UJ(i,i) = AUJ_UJ(i,i) - param.pf.lambda_o;
                    AUJ_OJ(i,i) = AUJ_OJ(i,i) + param.pf.lambda_o;
                elseif tremble ~= 0
                    t = tremble/exp(solutionSP.UJcValue(afi,ami,wmi) - solutionSP.OJcValue(afi,ami,wmi));
                    AUJ_UJ(i,i) = AUJ_UJ(i,i) - t;
                    AUJ_OJ(i,i) = AUJ_OJ(i,i) + t;
                end

                % UJ^c -> JJ^c (female job arrival over w_f)
                for wfi_new = 1:gridP.wf.n
                    jJJ = (wmi-1)*nam*na*gridP.wf.n + (ami-1)*na*gridP.wf.n + (wfi_new-1)*na + afi; % (af,wf,am,wm)
                    if solutionSP.UJcValue(afi,ami,wmi) < solutionSP.JJcValue(afi,wfi_new,ami,wmi)
                        AUJ_UJ(i,i)   = AUJ_UJ(i,i)   - param.pf.lambda_u * param.f(wfi_new);
                        AUJ_JJ(i,jJJ) = AUJ_JJ(i,jJJ) + param.pf.lambda_u * param.f(wfi_new);
                    elseif tremble ~= 0
                        t = tremble/exp(solutionSP.UJcValue(afi,ami,wmi) - solutionSP.JJcValue(afi,wfi_new,ami,wmi));
                        AUJ_UJ(i,i)   = AUJ_UJ(i,i)   - t * param.f(wfi_new);
                        AUJ_JJ(i,jJJ) = AUJ_JJ(i,jJJ) + t * param.f(wfi_new);
                    end
                end
            end
        end
    end

    AUJ_UJ = sparse(AUJ_UJ) + param.CUJ;
    AUJ_JJ = sparse(AUJ_JJ);
    AUJ_UU = param.CUJ_UU;
    AUJ_OJ = sparse(AUJ_OJ);
    AUJ_UO = sparse(AUJ_UO);
end

    function [AOJ_OJ, AOJ_OU, AOJ_OO, AOJ_JJ] = generateAOJc(solutionSP, param, gridP)
    AOJ_OJ = (zeros(param.NP/gridP.wf.n,param.NP/gridP.wf.n));    
    AOJ_OO = zeros(size(param.CUJ_UU));
    AOJ_JJ = (zeros(param.NP/gridP.wf.n,param.NP));
    AOJ_OU = zeros(size(param.CUJ_UU));

    for afi = 1:gridP.af.n
        for ami = 1:gridP.am.n
            for wmi = 1:gridP.wm.n
                i = (wmi-1)*gridP.af.n*gridP.am.n+(ami-1)*gridP.af.n+afi;                     

                % OJc -> OOc (male quits to OLF): unchanged
                i_new = (ami-1)*gridP.af.n+afi;        
                if solutionSP.OJcValue(afi, ami, wmi) < solutionSP.OOcValue(afi, ami)
                    AOJ_OJ(i,i)     = AOJ_OJ(i,i)     - param.pm.lambda_o;
                    AOJ_OO(i,i_new) = AOJ_OO(i,i_new) + param.pm.lambda_o;
                elseif tremble~=0
                    t = tremble/exp( solutionSP.OJcValue(afi, ami, wmi) - solutionSP.OOcValue(afi, ami) );
                    AOJ_OJ(i,i)     = AOJ_OJ(i,i)     - t;
                    AOJ_OO(i,i_new) = AOJ_OO(i,i_new) + t;
                end

                % NEW: OJc -> JJc via female OLF -> E (over wfi)
                for wfi_new = 1:gridP.wf.n
                    i_new_JJ = (wmi-1)*gridP.am.n*gridP.af.n*gridP.wf.n + ...
                               (ami-1)*gridP.af.n*gridP.wf.n + (wfi_new-1)*gridP.af.n + afi;
                    if solutionSP.OJcValue(afi, ami, wmi) < solutionSP.JJcValue(afi, wfi_new, ami, wmi)
                        AOJ_OJ(i,i)        = AOJ_OJ(i,i)        - param.pf.lambda_oE * param.f(wfi_new);
                        AOJ_JJ(i,i_new_JJ) = AOJ_JJ(i,i_new_JJ) + param.pf.lambda_oE * param.f(wfi_new);
                    elseif tremble~=0
                        t = tremble/exp( solutionSP.OJcValue(afi, ami, wmi) - solutionSP.JJcValue(afi, wfi_new, ami, wmi) );
                        AOJ_OJ(i,i)        = AOJ_OJ(i,i)        - t * param.f(wfi_new);
                        AOJ_JJ(i,i_new_JJ) = AOJ_JJ(i,i_new_JJ) + t * param.f(wfi_new);
                    end
                end
            end
        end        
    end 
    AOJ_OJ = sparse(AOJ_OJ)+ param.CUJ;
    AOJ_OU = param.CUJ_UU;
    AOJ_OO = sparse(AOJ_OO);
    AOJ_JJ = sparse(AOJ_JJ);
end

function [AUU_UU, AUU_UJ, AUU_JU, AUU_UO, AUU_OU] = generateAUUc(solutionSP, param, gridP)
    % UU state: (af,am)
    AUU_UU = zeros(gridP.af.n*gridP.am.n, gridP.af.n*gridP.am.n);
    AUU_UJ = zeros(gridP.af.n*gridP.am.n, param.NP/gridP.wf.n); % -> UJ (w_f)
    AUU_JU = zeros(gridP.af.n*gridP.am.n, param.NP/gridP.wm.n); % -> JU (w_m)
    AUU_UO = zeros(size(AUU_UU)); % -> UO
    AUU_OU = zeros(size(AUU_UU)); % -> OU

    for afi = 1:gridP.af.n
        for ami = 1:gridP.am.n
            i = (ami-1)*gridP.af.n + afi;

            % UU^c -> UO^c (male U->O)
            if solutionSP.UUcValue(afi,ami) < solutionSP.UOcValue(afi,ami)
                AUU_UU(i,i) = AUU_UU(i,i) - param.pm.lambda_o;
                AUU_UO(i,i) = AUU_UO(i,i) + param.pm.lambda_o;
            elseif tremble ~= 0
                t = tremble/exp(solutionSP.UUcValue(afi,ami) - solutionSP.UOcValue(afi,ami));
                AUU_UU(i,i) = AUU_UU(i,i) - t;
                AUU_UO(i,i) = AUU_UO(i,i) + t;
            end

            % UU^c -> OU^c (female U->O)
            if solutionSP.UUcValue(afi,ami) < solutionSP.OUcValue(afi,ami)
                AUU_UU(i,i) = AUU_UU(i,i) - param.pf.lambda_o;
                AUU_OU(i,i) = AUU_OU(i,i) + param.pf.lambda_o;
            elseif tremble ~= 0
                t = tremble/exp(solutionSP.UUcValue(afi,ami) - solutionSP.OUcValue(afi,ami));
                AUU_UU(i,i) = AUU_UU(i,i) - t;
                AUU_OU(i,i) = AUU_OU(i,i) + t;
            end

            % UU^c -> JU^c (female job arrival to JU over w_f)
            for wfi_new = 1:gridP.wf.n
                jJU = (ami-1)*gridP.af.n*gridP.wf.n + (wfi_new-1)*gridP.af.n + afi; % (af,wf,am)
                if solutionSP.UUcValue(afi,ami) < solutionSP.JUcValue(afi,wfi_new,ami)
                    AUU_UU(i,i) = AUU_UU(i,i) - param.pf.lambda_u * param.f(wfi_new);
                    AUU_JU(i,jJU) = AUU_JU(i,jJU) + param.pf.lambda_u * param.f(wfi_new);
                elseif tremble ~= 0
                    t = tremble/exp(solutionSP.UUcValue(afi,ami) - solutionSP.JUcValue(afi,wfi_new,ami));
                    AUU_UU(i,i) = AUU_UU(i,i) - t * param.f(wfi_new);
                    AUU_JU(i,jJU) = AUU_JU(i,jJU) + t * param.f(wfi_new);
                end
            end

            % UU^c -> UJ^c (male job arrival, over w_m)
            for wmi_new = 1:gridP.wm.n
                jUJ = (wmi_new-1)*gridP.af.n*gridP.am.n + (ami-1)*gridP.af.n + afi; % (af,am,wm)
                if solutionSP.UUcValue(afi,ami) < solutionSP.UJcValue(afi,ami,wmi_new)
                    AUU_UU(i,i) = AUU_UU(i,i) - param.pm.lambda_u * param.f(wmi_new);
                    AUU_UJ(i,jUJ) = AUU_UJ(i,jUJ) + param.pm.lambda_u * param.f(wmi_new);
                elseif tremble ~= 0
                    t = tremble/exp(solutionSP.UUcValue(afi,ami) - solutionSP.UJcValue(afi,ami,wmi_new));
                    AUU_UU(i,i) = AUU_UU(i,i) - t * param.f(wmi_new);
                    AUU_UJ(i,jUJ) = AUU_UJ(i,jUJ) + t * param.f(wmi_new);
                end
            end
        end
    end
    AUU_UU = sparse(AUU_UU) + param.CUU;
    AUU_UJ = sparse(AUU_UJ);
    AUU_JU = sparse(AUU_JU);
    AUU_UO = sparse(AUU_UO);
    AUU_OU = sparse(AUU_OU);
end


function [AOU_OU, AOU_OJ, AOU_OO, AOU_JU] = generateAOUc(solutionSP, param, gridP)
    AOU_OU = (zeros(gridP.af.n*gridP.am.n,gridP.af.n*gridP.am.n));
    AOU_OJ = (zeros(gridP.af.n*gridP.am.n,param.NP/gridP.wm.n));    
    AOU_OO = zeros(size(AOU_OU));
    AOU_JU = (zeros(gridP.af.n*gridP.am.n,param.NP/gridP.wm.n));

    for afi = 1:gridP.af.n
        for ami = 1:gridP.am.n
            i = (ami-1)*gridP.af.n+afi;   

            % OUc -> OOc (male U->O): unchanged
            if solutionSP.OUcValue(afi, ami) < solutionSP.OOcValue(afi, ami)
                AOU_OU(i,i) = AOU_OU(i,i) - param.pm.lambda_o;
                AOU_OO(i,i) = AOU_OO(i,i) + param.pm.lambda_o;
            elseif tremble~=0
                t = tremble/exp( solutionSP.OUcValue(afi, ami) - solutionSP.OOcValue(afi, ami) );
                AOU_OU(i,i) = AOU_OU(i,i) - t;
                AOU_OO(i,i) = AOU_OO(i,i) + t;
            end
            
            % NEW: OUc -> JUc via female OLF -> E (over wfi)
            for wfi_new = 1:gridP.wf.n
                i_new_JU = (ami-1)*gridP.af.n*gridP.wf.n + (wfi_new-1)*gridP.af.n + afi;
                if solutionSP.OUcValue(afi, ami) < solutionSP.JUcValue(afi, wfi_new, ami)
                    AOU_OU(i,i)      = AOU_OU(i,i)      - param.pf.lambda_oE * param.f(wfi_new);
                    AOU_JU(i,i_new_JU) = AOU_JU(i,i_new_JU) + param.pf.lambda_oE * param.f(wfi_new);
                elseif tremble~=0
                    t = tremble/exp( solutionSP.OUcValue(afi, ami) - solutionSP.JUcValue(afi, wfi_new, ami) );
                    AOU_OU(i,i)      = AOU_OU(i,i)      - t * param.f(wfi_new);
                    AOU_JU(i,i_new_JU) = AOU_JU(i,i_new_JU) + t * param.f(wfi_new);
                end
            end

            % male U->J (OUc -> OJc): unchanged
            for wmi_new = 1:gridP.wm.n
                i_new_OJ = (wmi_new-1)*gridP.af.n*gridP.am.n + (ami-1)*gridP.af.n + afi;
                if true&&(solutionSP.OUcValue(afi,ami)<solutionSP.OJcValue(afi,ami,wmi_new))
                    AOU_OU(i,i) = AOU_OU(i,i)-param.pm.lambda_u*param.f(wmi_new);
                    AOU_OJ(i,i_new_OJ) = AOU_OJ(i,i_new_OJ)+param.pm.lambda_u*param.f(wmi_new);
                elseif tremble~=0
                    t = tremble/exp(solutionSP.OUcValue(afi,ami)-solutionSP.OJcValue(afi,ami,wmi_new));
                    AOU_OU(i,i) = AOU_OU(i,i)-t*param.f(wmi_new);
                    AOU_OJ(i,i_new_OJ) = AOU_OJ(i,i_new_OJ)+t*param.f(wmi_new);
                end
            end
        end        
    end 
    AOU_OU = sparse(AOU_OU) + param.CUU;
    AOU_OJ = sparse(AOU_OJ);
    AOU_OO = sparse(AOU_OO);
    AOU_JU = sparse(AOU_JU);
end

function [AUO_UO, AUO_JO, AUO_OO, AUO_UJ] = generateAUOc(solutionSP, param, gridP)
    AUO_UO = (zeros(gridP.af.n*gridP.am.n,gridP.af.n*gridP.am.n));    
    AUO_JO = (zeros(gridP.af.n*gridP.am.n,param.NP/gridP.wm.n));
    AUO_OO = zeros(size(AUO_UO));
    AUO_UJ = (zeros(gridP.af.n*gridP.am.n,param.NP/gridP.wf.n));

    for afi = 1:gridP.af.n
        for ami = 1:gridP.am.n
            i = (ami-1)*gridP.af.n+afi;

            % UOc -> OOc (female U->O): unchanged
            if solutionSP.UOcValue(afi, ami) < solutionSP.OOcValue(afi, ami)
                AUO_UO(i,i) = AUO_UO(i,i) - param.pf.lambda_o;
                AUO_OO(i,i) = AUO_OO(i,i) + param.pf.lambda_o;
            elseif tremble~=0
                t = tremble/exp( solutionSP.UOcValue(afi, ami) - solutionSP.OOcValue(afi, ami) );
                AUO_UO(i,i) = AUO_UO(i,i) - t;
                AUO_OO(i,i) = AUO_OO(i,i) + t;
            end
            
            % NEW: UOc -> UJc via male OLF -> E (over wmi)
            for wmi_new = 1:gridP.wm.n
                i_new_UJ = (wmi_new-1)*gridP.af.n*gridP.am.n + (ami-1)*gridP.af.n + afi;
                if solutionSP.UOcValue(afi, ami) < solutionSP.UJcValue(afi, ami, wmi_new)
                    AUO_UO(i,i)       = AUO_UO(i,i)       - param.pm.lambda_oE * param.f(wmi_new);
                    AUO_UJ(i,i_new_UJ) = AUO_UJ(i,i_new_UJ) + param.pm.lambda_oE * param.f(wmi_new);
                elseif tremble~=0
                    t = tremble/exp( solutionSP.UOcValue(afi, ami) - solutionSP.UJcValue(afi, ami, wmi_new) );
                    AUO_UO(i,i)       = AUO_UO(i,i)       - t * param.f(wmi_new);
                    AUO_UJ(i,i_new_UJ) = AUO_UJ(i,i_new_UJ) + t * param.f(wmi_new);
                end
            end

            % female U->J (UOc -> JOc): unchanged
            for wfi_new = 1:gridP.wf.n
                i_new_JO = (ami-1)*gridP.af.n*gridP.wf.n+(wfi_new-1)*gridP.af.n+afi;
                if true&&(solutionSP.UOcValue(afi,ami)<solutionSP.JOcValue(afi,wfi_new,ami))
                    AUO_UO(i,i) = AUO_UO(i,i)-param.pf.lambda_u*param.f(wfi_new);
                    AUO_JO(i,i_new_JO) = AUO_JO(i,i_new_JO)+param.pf.lambda_u*param.f(wfi_new);
                elseif tremble~=0
                    t = tremble/exp(solutionSP.UOcValue(afi,ami)-solutionSP.JOcValue(afi,wfi_new,ami));
                    AUO_UO(i,i) = AUO_UO(i,i)-t*param.f(wfi_new);
                    AUO_JO(i,i_new_JO) = AUO_JO(i,i_new_JO)+t*param.f(wfi_new);
                end
            end
        end        
    end 
    AUO_UO = sparse(AUO_UO) + param.CUU;
    AUO_JO = sparse(AUO_JO);
    AUO_OO = sparse(AUO_OO);
    AUO_UJ = sparse(AUO_UJ);
end

function [AOO_OO, AOO_OJ, AOO_JO] = generateAOOc(solutionSP, param, gridP)
    AOO_OO = (zeros(gridP.af.n*gridP.am.n,gridP.af.n*gridP.am.n));        
    AOO_OJ = (zeros(gridP.af.n*gridP.am.n, param.NP/gridP.wf.n)); % male to E -> OJc
    AOO_JO = (zeros(gridP.af.n*gridP.am.n, param.NP/gridP.wm.n)); % female to E -> JOc

    for afi = 1:gridP.af.n
        for ami = 1:gridP.am.n
            i = (ami-1)*gridP.af.n+afi; 

            % NEW: OOc -> JOc via female OLF -> E (over wfi)
            for wfi_new = 1:gridP.wf.n
                i_new_JO = (ami-1)*gridP.af.n*gridP.wf.n + (wfi_new-1)*gridP.af.n + afi;
                if solutionSP.OOcValue(afi, ami) < solutionSP.JOcValue(afi, wfi_new, ami)
                    AOO_OO(i,i)      = AOO_OO(i,i)      - param.pf.lambda_oE * param.f(wfi_new);
                    AOO_JO(i,i_new_JO) = AOO_JO(i,i_new_JO) + param.pf.lambda_oE * param.f(wfi_new);
                elseif tremble~=0
                    t = tremble/exp( solutionSP.OOcValue(afi, ami) - solutionSP.JOcValue(afi, wfi_new, ami) );
                    AOO_OO(i,i)      = AOO_OO(i,i)      - t * param.f(wfi_new);
                    AOO_JO(i,i_new_JO) = AOO_JO(i,i_new_JO) + t * param.f(wfi_new);
                end
            end
            
            % NEW: OOc -> OJc via male OLF -> E (over wmi)
            for wmi_new = 1:gridP.wm.n
                i_new_OJ = (wmi_new-1)*gridP.af.n*gridP.am.n + (ami-1)*gridP.af.n + afi;
                if solutionSP.OOcValue(afi, ami) < solutionSP.OJcValue(afi, ami, wmi_new)
                    AOO_OO(i,i)      = AOO_OO(i,i)      - param.pm.lambda_oE * param.f(wmi_new);
                    AOO_OJ(i,i_new_OJ) = AOO_OJ(i,i_new_OJ) + param.pm.lambda_oE * param.f(wmi_new);
                elseif tremble~=0
                    t = tremble/exp( solutionSP.OOcValue(afi, ami) - solutionSP.OJcValue(afi, ami, wmi_new) );
                    AOO_OO(i,i)      = AOO_OO(i,i)      - t * param.f(wmi_new);
                    AOO_OJ(i,i_new_OJ) = AOO_OJ(i,i_new_OJ) + t * param.f(wmi_new);
                end
            end
        end        
    end 
    AOO_OO = sparse(AOO_OO) + param.CUU;    
    AOO_OJ = sparse(AOO_OJ);
    AOO_JO = sparse(AOO_JO);
end


% =========================================================================
% Instantaneous utilities (paper, Section 2.1)
% =========================================================================
% Singles: c(W,k) = D_k(w e^a) - k kappa(a),  c(U,k) = D_k(b(a)) - k kappa(a),
%          c(O,k) = h^s_{g,k} - k kappa(a).
% Couples: c = [D_k(y_f + y_m) + H_f + H_m - k kappa((a_f + a_m)/2)] / 2 per
%          spouse, with y = w e^a (W), b(a) (U), 0 (O) and H = h^c_{g,k} (O).
% D_0 = param.nokid_tax, D_1 = param.onekid_tax; home production is not
% taxed. Women's wage grids already include theta_w (Build*Grid).

function u = instUtilityOl(gridS,pg)
    u = util(param, pg.hs0 .* ones(gridS.a.n,1));
end
function u = instUtilityUl(gridS,pg) %#ok<INUSD>
    u = util(param, param.nokid_tax(param.b(gridS.a.a),gridS.a.a,param));
end
function u = instUtilityJl(gridS,pg) %#ok<INUSD>
    u = util(param, param.nokid_tax(gridS.w.ww.*exp(gridS.a.aa),gridS.a.a,param));
end
function u = instUtilityOc(gridS,pg)
    u = util(param, pg.hs1 - param.c(gridS.a.a));
end
function u = instUtilityUc(gridS,pg) %#ok<INUSD>
    u = util(param, param.onekid_tax(param.b(gridS.a.a),gridS.a.a,param) - param.c(gridS.a.a));
end
function u = instUtilityJc(gridS,pg) %#ok<INUSD>
    u = util(param, param.onekid_tax(gridS.w.ww.*exp(gridS.a.aa),gridS.a.a,param) - param.c(gridS.a.aa));
end

function u = coupleUtility(y, H, k)
    % Per-spouse utility on the [af x wf x am x wm] grid: y = pooled taxable
    % income, H = home production of spouses out of the labor force,
    % k = child indicator. One household child cost at average productivity.
    amid = (gridP.af.aa + gridP.am.aa)./2;
    y    = y + zeros(size(amid));
    if k
        cons = (param.onekid_tax(y, amid, param) + H - param.c(amid)) ./ 2;
    else
        cons = (param.nokid_tax(y, amid, param) + H) ./ 2;
    end
    u = util(param, cons);
end

% =========================
% Couple utilities - no kid
% =========================

function u = instUtilityJJl()
    Yf = gridP.wf.ww .* exp(gridP.af.aa);
    Ym = gridP.wm.ww .* exp(gridP.am.aa);
    u  = coupleUtility(Yf + Ym, 0, 0);                         % [af x wf x am x wm]
end
function u = instUtilityJUl()
    Yf = gridP.wf.ww .* exp(gridP.af.aa);
    u4 = coupleUtility(Yf + param.b(gridP.am.aa), 0, 0);
    u  = squeeze(u4(:,:,:,1));                                 % [af x wf x am]
end
function u = instUtilityUJl()
    Ym = gridP.wm.ww .* exp(gridP.am.aa);
    u4 = coupleUtility(param.b(gridP.af.aa) + Ym, 0, 0);
    u  = squeeze(u4(:,1,:,:));                                 % [af x am x wm]
end
function u = instUtilityUUl()
    u4 = coupleUtility(param.b(gridP.af.aa) + param.b(gridP.am.aa), 0, 0);
    u  = squeeze(u4(:,1,:,1));                                 % [af x am]
end
function u = instUtilityJOl()
    Yf = gridP.wf.ww .* exp(gridP.af.aa);
    u4 = coupleUtility(Yf, param.pm.hc0, 0);
    u  = squeeze(u4(:,:,:,1));                                 % [af x wf x am]
end
function u = instUtilityOJl()
    Ym = gridP.wm.ww .* exp(gridP.am.aa);
    u4 = coupleUtility(Ym, param.pf.hc0, 0);
    u  = squeeze(u4(:,1,:,:));                                 % [af x am x wm]
end
function u = instUtilityUOl()
    u4 = coupleUtility(param.b(gridP.af.aa), param.pm.hc0, 0);
    u  = squeeze(u4(:,1,:,1));                                 % [af x am]
end
function u = instUtilityOUl()
    u4 = coupleUtility(param.b(gridP.am.aa), param.pf.hc0, 0);
    u  = squeeze(u4(:,1,:,1));                                 % [af x am]
end
function u = instUtilityOOl()
    u4 = coupleUtility(0, param.pf.hc0 + param.pm.hc0, 0);
    u  = squeeze(u4(:,1,:,1));                                 % [af x am]
end

% ========================
% Couple utilities - 1 kid
% ========================

function u = instUtilityJJc()
    Yf = gridP.wf.ww .* exp(gridP.af.aa);
    Ym = gridP.wm.ww .* exp(gridP.am.aa);
    u  = coupleUtility(Yf + Ym, 0, 1);                         % [af x wf x am x wm]
end
function u = instUtilityJUc()
    Yf = gridP.wf.ww .* exp(gridP.af.aa);
    u4 = coupleUtility(Yf + param.b(gridP.am.aa), 0, 1);
    u  = squeeze(u4(:,:,:,1));                                 % [af x wf x am]
end
function u = instUtilityUJc()
    Ym = gridP.wm.ww .* exp(gridP.am.aa);
    u4 = coupleUtility(param.b(gridP.af.aa) + Ym, 0, 1);
    u  = squeeze(u4(:,1,:,:));                                 % [af x am x wm]
end
function u = instUtilityUUc()
    u4 = coupleUtility(param.b(gridP.af.aa) + param.b(gridP.am.aa), 0, 1);
    u  = squeeze(u4(:,1,:,1));                                 % [af x am]
end
function u = instUtilityJOc()
    Yf = gridP.wf.ww .* exp(gridP.af.aa);
    u4 = coupleUtility(Yf, param.pm.hc1, 1);
    u  = squeeze(u4(:,:,:,1));                                 % [af x wf x am]
end
function u = instUtilityOJc()
    Ym = gridP.wm.ww .* exp(gridP.am.aa);
    u4 = coupleUtility(Ym, param.pf.hc1, 1);
    u  = squeeze(u4(:,1,:,:));                                 % [af x am x wm]
end
function u = instUtilityUOc()
    u4 = coupleUtility(param.b(gridP.af.aa), param.pm.hc1, 1);
    u  = squeeze(u4(:,1,:,1));                                 % [af x am]
end
function u = instUtilityOUc()
    u4 = coupleUtility(param.b(gridP.am.aa), param.pf.hc1, 1);
    u  = squeeze(u4(:,1,:,1));                                 % [af x am]
end
function u = instUtilityOOc()
    u4 = coupleUtility(0, param.pf.hc1 + param.pm.hc1, 1);
    u  = squeeze(u4(:,1,:,1));                                 % [af x am]
end

end