function mapPoints = getMapPoints(dataMap, startDate, endDate)
%  function mapPoints = getMapPoints(dataMap, startDate, endDate)
% return an array of map points which have data within the listed period
% from startDate to endDate (Matlab date system). 
mapPoints = [];
for i = 1:length(dataMap);
   aPoint = dataMap(i);
   if aPoint.endDate < startDate || aPoint.startDate > endDate
       continue;
   end
   if isempty(mapPoints)
       mapPoints = aPoint;
   else
       mapPoints = [mapPoints aPoint];
   end
end