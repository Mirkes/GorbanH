function saveFigures(fName, drive, pos)
    set(gcf, 'PaperPositionMode', 'auto');
    set(gcf,'color','w');
    if nargin > 2
        set(gcf, 'pos', pos);
    end
    print(gcf, drive, '-noui', '-loose', fName);
end