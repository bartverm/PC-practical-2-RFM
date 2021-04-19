function plot_zb_sym(x,z_b_sym, dt_sym, t_end, So)
t_unit='s';
if nargin < 5
    So=0;
    if nargin < 4
        t_end=5;
        if nargin < 3
            dt_sym=1;
            t_unit='';
        end
    end
end
nt=size(z_b_sym,2);
dt=t_end/nt;
z_b_sym=z_b_sym+((x-max(x))*So)';
hf=gcf;
hbed=plot(x,z_b_sym(:,1),'k','linewidth',2);
ht=title(['t = ',num2str(dt_sym,3),' ',t_unit]);
set(gca,'ylim',[min(z_b_sym(:)) max(z_b_sym(:))],'ylimmode','manual')
xlabel('x (m)')
ylabel('z (m)')
ct=2;
while isvalid(hf) && ct <=nt
    hbed.YData=z_b_sym(:,ct);
    ht.String=['t = ',num2str(ct*dt_sym,3),' ',t_unit];
    pause(dt)
    ct=ct+1;
end
end