
% Export DIRAC results to Excel

load BeeDIRACComps_KEGG_100213
sheetnames = fieldnames(results_tables);

for i = 1:numel(comparison)
    filename = [comparison{i},'_KEGG.xslx'];
    for j = 1:numel(sheetnames)
        sheetname = sheetnames{j};
        table = getfield(results_tables,sheetname);
        xlswrite(filename,table,sheetname);
    end
end
    