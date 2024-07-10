function [] = plot_globe()

    load('topo.mat','topo');
    
    R_E = 6378e3;
    [px,py,pz] = sphere;
    pxx = px * R_E;
    pyy = py * R_E;
    pzz = pz * R_E;

    sEarth = surface(pyy, -pxx , pzz);
    sEarth.FaceColor = 'texturemap';
    sEarth.EdgeColor = '[0.35,0.35,0.35]';
    sEarth.LineStyle = ":";
    sEarth.CData = topo;
   
    daspect([1 1 1])
    view([103 25])

    axis tight

    % Rotate the map
    % rotate(sEarth, [0 0 1], 120.9605); % Rotating 30 degrees around the z-axis

end