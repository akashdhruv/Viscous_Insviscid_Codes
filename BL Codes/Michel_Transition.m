function [x_tp, z_tp, s_tp] = Michel_Transition(Re_x, Re_d2, x, z, s)

MT_t = 1.174*(1+(22400./Re_x)).*(Re_x.^0.46);
MT_t(1) = 0;

Input1 = MT_t-Re_d2;

%         figure;
%         plot(s,Input1,'-k')
%         xlabel('Arc Length')
%         ylabel('Result')
        
input = pchip(s,Input1);

s_tp = fzero(@fn1,0.5);                                                    % Arc Length Location of the Transition Point

    function u = fn1(s)
        u = ppval(input,s);
    end

% % % j = length(MT_t);
% % %    
% % % ind = zeros(length(Re_x),1);
% % % 
% % % for k = 1:j;
% % %     if gt(Re_d2(k),MT_t(k));
% % %         ind(k) = 1;
% % %     else
% % %         ind(k) = 0;
% % %     end
% % % end     
% % % 
% % % ffz = Re_d2-MT_t;                                                          % Function to find zero - The Transition Point
% % % 
% % % s_tp = ppval(pchip(ffz(2:end),s(2:end)),0);                                % Transition Point Arc Length
% % %     
% % % Re_t_T = ppval(pchip(s,Re_d2),s_tp);                                       % Momentum Thickness Dependent Re at the Transition Point

%     figure;
%     plot(s,Re_d2,'-k',s,MT_t,'-g')
%     title('Michels Transition Point Identification')
%     xlabel('Normalized Arc Length')
%     ylabel('Momentum Thickness Dependant Reynolds Number, Re_\theta') 

if isnan(s_tp)
    s_tp = s(end);
end

Re_d2_TP = ppval(pchip(s,Re_d2),s_tp);
Re_x_TP  = ppval(pchip(s,Re_x), s_tp);

x_tp = ppval(pchip(s,x),s_tp);
z_tp = ppval(pchip(s,z),s_tp);
end