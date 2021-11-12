function plot_labels_diff(ypos, max_temp, S)

delta = 0.125;

axis([min(S) max(S) -1 max_temp])
text((S(2)-S(1))/2 + S(1)-delta,ypos,'$\mathcal{N}_1$','interpreter','latex','fontsize',20)
text((S(3)-S(2))/2 + S(2)-delta,ypos,'$\mathcal{N}_2$','interpreter','latex','fontsize',20)
text((S(4)-S(3))/2 + S(3)-delta,ypos,'$\mathcal{N}_3$','interpreter','latex','fontsize',20)

end