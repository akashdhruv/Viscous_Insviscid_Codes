function [du_ds] = vel_grad(s,ue)

% for i = 1:(length(sj)-1)
%     aem(i) = (uej(i+1)-uej(i))/(sj(i+1)-sj(i));                            % Velocity Gradient at Panel Midpoints
%     sm(i) = sj(i)+(sj(i+1)-sj(i))/2;                                       % Arc Length Location of Velocity at Midpoints
% end
% 
% %         figure;
% %         plot(sm,aem,'-k')
% %         title('Velocity Gradient at Panel Midpoints')
% %         xlabel('Arc Length')
% %         ylabel('Velocity Gradient')
%         
% aej = ppval(pchip(sm,aem),sj);                                      
    

for i = 1:length(s)
        if i == 1
            ds(i) = s(2)-s(1);
            du = ue(2)-ue(1);
            du_ds(:,i) = du/ds(i);                                            % Change in Velocity wrt to Arc Length
        elseif i == length(s)
            ds(i) = s(end-1)-s(end);
            du = ue(end-1)-ue(end);
            du_ds(:,i) = du/ds(i);
        else
            ds(i) = s(i)-s(i-1);
            du = ue(i)-ue(i-1);
            du_ds(:,i) = du/ds(i);
        end
end
end

