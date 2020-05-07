function commas
     xl = get(gca,'XTickLabel');
     new_xl = strrep(xl(:),'.',',');
     set(gca,'XTickLabel',new_xl);

     yl = get(gca,'YTickLabel');
     new_yl = strrep(yl(:),'.',',');
     set(gca,'YTickLabel',new_yl);
end

