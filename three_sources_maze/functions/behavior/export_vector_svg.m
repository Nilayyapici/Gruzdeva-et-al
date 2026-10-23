function export_vector_svg(fig, filename)
    if nargin < 1 || isempty(fig)
        fig = gcf;
    end
    if nargin < 2
        filename = 'figure_export';
    end
    
    drawnow;
    all_patches = findobj(fig, 'Type', 'patch');
    
    if ~isempty(all_patches)
        original_alphas = arrayfun(@(x) x.FaceAlpha, all_patches);
        arrayfun(@(x) set(x, 'FaceAlpha', 1), all_patches);
        drawnow;
    end
    
    print(fig, filename, '-dsvg', '-painters');
    fprintf('Exported vector SVG: %s.svg\n', filename);
    
    if ~isempty(all_patches)
        for i = 1:length(all_patches)
            all_patches(i).FaceAlpha = original_alphas(i);
        end
        drawnow;
    end
end