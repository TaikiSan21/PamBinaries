function millis = dateNumToMillis(datenum)

millis = (datenum-719529)*86400000;
% datenum = double(millis)/86400000.0+719529;