function Y = save_pdf_without_whitespace3(name);
set(gca,'FontSize',20);
fig = gcf;
fig.PaperPositionMode = 'auto';
print(fig,name,'-dpdf');
end
