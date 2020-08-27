function[map] = nova(N)
% N - number of colours in the map
vec = [100:-10:0];
% Nova
hex = ['#4706ac';'#450cb3';'#4212bb';'#3a1cca';'#a900ae';'#c600a0';'#ed007a';'#fd3c54';'#fb6b31';'#eb940e';'#FF9B00']; 
raw = sscanf(hex','#%2x%2x%2x',[3,size(hex,1)]).' / 255;
map = interp1(vec,raw,linspace(100,0,N),'pchip');
% %% Check map in RGB
% figure;
% rgbplot(map);
% colormap(map);
% colorbar('Ticks',[]);
% figure;
% colormap(map);
% for iCol = 1:size(map,1)
% scatter3(map(iCol,1), map(iCol,2), map(iCol,3), 50, map(iCol,:), 'filled'); hold on;
% end
% xlabel('Red'); ylabel('Green'); zlabel('Blue');
% colorbar('Ticks',[]);
end

